#!/usr/bin/env python3
"""Count perfect matches the way the SUNi paper does.

Mighell, Toledano, and Lehner, PLOS ONE 2023, count a read only when the merged
sequence is an exact match to an expected NNK or NNS allele. Upstream of this
script, vsearch merges the pair and discards merges with expected error >= 0.5,
and cutadapt trims a 3' Illumina adapter when one is present.

Reads FASTQ on stdin. The alleles JSON is --meta. The summary path is the
positional argument.

    python count_paper_style.py --meta designed_alleles.json summary.json
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

from analyze_bzip_library import CODON_AA

ROOT = Path(__file__).resolve().parent
META = ROOT / "bzip_mutation_profile" / "designed_alleles.json"

DEGENERATE = {"N": "ACGT", "K": "GT", "S": "GC"}


def expand_codon(pattern: str) -> set[str]:
    options = {""}
    for base in pattern:
        alphabet = DEGENERATE.get(base, base)
        options = {prefix + nt for prefix in options for nt in alphabet}
    return options


def percentile(sorted_vals: list[float], percent: float) -> float:
    """Linear percentile, matching numpy.percentile default."""
    if not sorted_vals:
        raise ValueError("empty")
    k = (len(sorted_vals) - 1) * percent / 100.0
    lo = math.floor(k)
    hi = math.ceil(k)
    if lo == hi:
        return sorted_vals[lo]
    return sorted_vals[lo] * (hi - k) + sorted_vals[hi] * (k - lo)


def load_reference(meta_path: Path) -> tuple[str, list[int | None], dict[int, int], dict[int, set[str]], dict[int, str]]:
    meta = json.loads(meta_path.read_text())
    wt = meta["amplicon"]
    amp_len = int(meta["amplicon_length"])
    if len(wt) != amp_len:
        raise SystemExit(f"amplicon length {len(wt)} != {amp_len}")
    codon_at: list[int | None] = [None] * amp_len
    offset_of: dict[int, int] = {}
    allowed: dict[int, set[str]] = {}
    wt_codon: dict[int, str] = {}
    for key, info in meta["codons"].items():
        if not info.get("in_amplicon"):
            continue
        pos = int(key)
        offset = int(info["amp_offset"])
        pattern = info["degeneracy"]
        codons = expand_codon(pattern)
        observed = wt[offset : offset + 3]
        if observed != info["wt_codon"]:
            raise SystemExit(f"codon {pos} WT {observed} != {info['wt_codon']}")
        allowed[pos] = codons
        wt_codon[pos] = observed
        offset_of[pos] = offset
        for i in range(3):
            codon_at[offset + i] = pos
    return wt, codon_at, offset_of, allowed, wt_codon


def classify(
    seq: str,
    wt: str,
    codon_at: list[int | None],
    offset_of: dict[int, int],
    allowed: dict[int, set[str]],
    wt_codon: dict[int, str],
) -> tuple[str, int | None, str | None]:
    """Return (class, codon position, triplet). class is wt, programmed, or other."""
    if len(seq) != len(wt):
        return "other", None, None
    if seq == wt:
        return "wt", None, None
    changed: dict[int, int] = {}
    for i, (base, ref) in enumerate(zip(seq, wt)):
        if base == ref:
            continue
        pos = codon_at[i]
        if pos is None:
            return "other", None, None
        changed[pos] = changed.get(pos, 0) + 1
    if len(changed) != 1:
        return "other", None, None
    pos = next(iter(changed))
    triplet = seq[offset_of[pos] : offset_of[pos] + 3]
    if triplet == wt_codon[pos] or triplet not in allowed[pos]:
        return "other", None, None
    return "programmed", pos, triplet


def main() -> None:
    parser = argparse.ArgumentParser(prog="count_paper_style.py")
    parser.add_argument("--meta", type=Path, default=META)
    parser.add_argument("summary", type=Path)
    args = parser.parse_args()
    out_path = args.summary
    wt, codon_at, offset_of, allowed, wt_codon = load_reference(args.meta)
    positions = sorted(allowed)
    classes: Counter[str] = Counter()
    per_codon: Counter[int] = Counter()
    alleles: Counter[tuple[int, str]] = Counter()
    lengths: Counter[int] = Counter()
    n = 0
    header = sys.stdin.readline()
    while header:
        seq = sys.stdin.readline().strip().upper()
        plus = sys.stdin.readline()
        qual = sys.stdin.readline()
        if not qual:
            raise SystemExit("truncated FASTQ")
        if not header.startswith("@") or not plus.startswith("+"):
            raise SystemExit("FASTQ record out of frame")
        label, pos, triplet = classify(seq, wt, codon_at, offset_of, allowed, wt_codon)
        classes[label] += 1
        lengths[len(seq)] += 1
        if pos is not None and triplet is not None:
            per_codon[pos] += 1
            alleles[(pos, triplet)] += 1
        n += 1
        if n % 5_000_000 == 0:
            print(f"counted {n:,}", file=sys.stderr, flush=True)
        header = sys.stdin.readline()

    retained = n
    programmed = classes["programmed"]
    wildtype = classes["wt"]
    freqs = [per_codon[pos] / retained if retained else 0.0 for pos in positions]
    ordered = sorted(freqs)
    p10 = percentile(ordered, 10) if retained else 0.0
    p90 = percentile(ordered, 90) if retained else 0.0
    logdiff = math.log10(p90 / p10) if p10 > 0 and p90 > 0 else None
    programmed_fraction = programmed / retained if retained else 0.0
    efficiency = (
        programmed_fraction / (10 ** logdiff) if logdiff is not None else None
    )
    by_pos: dict[int, list[tuple[str, int]]] = defaultdict(list)
    for (pos, triplet), count in alleles.items():
        by_pos[pos].append((triplet, count))
    codon_rows = []
    for pos in positions:
        items = sorted(by_pos[pos], key=lambda item: -item[1])
        amino_acids = {CODON_AA[triplet] for triplet, _ in items}
        wt_aa = CODON_AA[wt_codon[pos]]
        codon_rows.append({
            "position": pos,
            "wt_codon": wt_codon[pos],
            "wt_aa": wt_aa,
            "count": per_codon[pos],
            "frequency": per_codon[pos] / retained if retained else 0.0,
            "non_wt_aa": len(amino_acids - {"*", wt_aa}),
            "stop_count": sum(count for triplet, count in items if CODON_AA[triplet] == "*"),
            "designed_seen": len(items) + (1 if wt_codon[pos] in allowed[pos] else 0),
            "top": [
                {"codon": triplet, "aa": CODON_AA[triplet], "count": count}
                for triplet, count in items[:8]
            ],
        })
    top_alleles = [
        {
            "position": pos,
            "codon": triplet,
            "aa": CODON_AA[triplet],
            "wt_aa": CODON_AA[wt_codon[pos]],
            "count": count,
        }
        for (pos, triplet), count in alleles.most_common(8)
    ]
    summary = {
        "method": (
            "vsearch --fastq_mergepairs --fastq_maxee 0.5, "
            "cutadapt -a AGATCGGAAGAGC --overlap 10, "
            "perfect match to an expected NNK/NNS sequence"
        ),
        "retained_reads": retained,
        "classes": dict(classes),
        "programmed_fraction": programmed_fraction,
        "wildtype_fraction": wildtype / retained if retained else 0.0,
        "other_fraction": classes["other"] / retained if retained else 0.0,
        "percentile_10": p10,
        "percentile_90": p90,
        "logdiff": logdiff,
        "screening_efficiency": efficiency,
        "per_codon": {str(pos): per_codon[pos] for pos in positions},
        "codons": codon_rows,
        "top_alleles": top_alleles,
        "alleles": [
            {
                "position": pos,
                "codon": triplet,
                "aa": CODON_AA[triplet],
                "count": count,
            }
            for (pos, triplet), count in alleles.most_common()
        ],
        "length_top": lengths.most_common(8),
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: summary[k] for k in (
        "retained_reads", "classes", "programmed_fraction", "wildtype_fraction",
        "other_fraction", "logdiff", "screening_efficiency", "percentile_10", "percentile_90",
    )}, indent=2))


if __name__ == "__main__":
    main()
