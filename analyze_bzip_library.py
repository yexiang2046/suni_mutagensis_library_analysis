#!/usr/bin/env python3
"""Designed-allele reference and bowtie2 tally for the JunB bZIP SUNI library.

Build the reference (wild-type amplicon plus every NNK/NNS allele):
    python analyze_bzip_library.py build

Map pairs with bowtie2 against the wild-type amplicon and tally:
    bowtie2 --end-to-end --very-sensitive -p <threads> \\
        -x bzip_mutation_profile/wt_amplicon \\
        -1 R1.fastq.gz -2 R2.fastq.gz \\
        | python analyze_bzip_library.py tally

Bowtie2 2.5.5 `--best` aborts with internal exception #1 on these pairs, and a
multi-FASTA of the one-codon alleles returns an arbitrary near-match (the
references differ by a single codon, so the FM-index cannot rank them). The
tally therefore aligns each pair end-to-end to the wild-type amplicon, rebuilds
the 211 bp read, and scores that sequence against every designed allele. An
exact call is the unique allele at edit distance 0. A tie is ambiguous. A
unique nearest allele that still has edits is an extra-edit.
"""

from __future__ import annotations

import csv
import json
import re
import sys
import zipfile
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "bzip_mutation_profile"
CDNA_DOCX = ROOT / "JunB cDNA.docx"
PRIMER_CSV = ROOT / "JunB SUNI Primers" / "primers_output_JunB_bZIP_codons267-330.csv"

ANCHOR = "ATGTGCACTAAAATGGAACAGCCCTTCTAC"
AMP_START = 1062
AMP_LEN = 211
READ_LEN = 151
# R2 reverse-complements onto this 0-based offset of the amplicon.
R2_OFFSET = 60
ANCHOR_5P = "CATCAACATGGAAGACCAAG"
DESIGN_FIRST = 267
DESIGN_LAST = 329
FLANK_FIRST = 260
FLANK_LAST = 266

DEGENERATE = {"N": "ACGT", "K": "GT", "S": "GC"}
CODON_AA = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY*")


def extract_cdna() -> str:
    xml = zipfile.ZipFile(CDNA_DOCX).read("word/document.xml").decode()
    text = re.sub(r"<[^>]+>", "", re.sub(r"</w:p>", "\n", xml))
    seqs = re.findall(r"[acgtACGT]{100,}", text)
    if not seqs:
        raise SystemExit(f"No cDNA sequence found in {CDNA_DOCX}")
    return max(seqs, key=len).upper()


def expand_codon(pattern: str) -> list[str]:
    options = [""]
    for base in pattern:
        alphabet = DEGENERATE.get(base, base)
        options = [prefix + nt for prefix in options for nt in alphabet]
    return options


def codon_genomic_start(anchor_index: int, codon: int) -> int:
    return anchor_index + (codon - 1) * 3


def amplicon_codon_offset(codon: int) -> int:
    """0-based offset of a codon inside the 211 bp amplicon."""
    return 1 + (codon - 260) * 3


def load_primers() -> list[dict]:
    primers = []
    with PRIMER_CSV.open() as handle:
        for row in csv.DictReader(handle):
            label = row["Primer Number"].strip()
            match = re.fullmatch(r"Primer\s+(\d+)", label)
            if not match:
                raise SystemExit(f"Unexpected primer label: {label}")
            sequence = row["Primer Sequence"].strip().upper()
            degenerates = re.findall(r"NNK|NNS", sequence)
            if len(degenerates) != 1:
                raise SystemExit(f"{label} has {degenerates}, expected one NNK/NNS codon")
            primers.append({
                "position": int(match.group(1)),
                "sequence": sequence,
                "degeneracy": degenerates[0],
            })
    return primers


def build() -> None:
    cdna = extract_cdna()
    anchor_index = cdna.find(ANCHOR)
    if anchor_index < 0:
        raise SystemExit("ORF anchor not found in the cDNA")
    if cdna[AMP_START:AMP_START + 20] != ANCHOR_5P:
        raise SystemExit(
            f"Amplicon start {AMP_START} is {cdna[AMP_START:AMP_START + 20]}, expected {ANCHOR_5P}"
        )
    amplicon = cdna[AMP_START:AMP_START + AMP_LEN]
    if len(amplicon) != AMP_LEN:
        raise SystemExit(f"Amplicon is {len(amplicon)} bp, expected {AMP_LEN}")

    # Codon 260 starts one base into the amplicon; codon 329 ends on the last base.
    if amplicon_codon_offset(260) != 1 or amplicon_codon_offset(329) + 3 != AMP_LEN:
        raise SystemExit("Codon frame does not fill the amplicon as codons 260-329")

    primers = load_primers()
    by_pos = {item["position"]: item for item in primers}
    expected = list(range(DESIGN_FIRST, DESIGN_LAST + 2))  # 267-330 inclusive
    if sorted(by_pos) != expected:
        raise SystemExit(f"Primer positions {sorted(by_pos)} != {expected[0]}-{expected[-1]}")

    OUT.mkdir(parents=True, exist_ok=True)
    alleles: dict[str, str] = {"WT": amplicon}
    annotations = {
        "WT": {"position": None, "codon": None, "aa": None, "degeneracy": None}
    }
    codons = {}
    skipped_wt = []
    skipped_outside = []

    for position in expected:
        primer = by_pos[position]
        genomic = codon_genomic_start(anchor_index, position)
        wt_codon = cdna[genomic:genomic + 3]
        if len(wt_codon) != 3 or set(wt_codon) - set("ACGT"):
            raise SystemExit(f"Codon {position} is not a DNA triplet: {wt_codon}")
        slot = primer["sequence"].find(primer["degeneracy"])
        concrete = primer["sequence"].replace(primer["degeneracy"], wt_codon, 1)
        placed = cdna.find(concrete)
        if placed < 0 or placed + slot != genomic:
            raise SystemExit(
                f"Primer {position} does not sit on codon {position} "
                f"(placed {placed}, slot {slot}, genomic {genomic})"
            )

        offset = genomic - AMP_START
        inside = 0 <= offset and offset + 3 <= AMP_LEN
        codons[str(position)] = {
            "wt_codon": wt_codon,
            "wt_aa": CODON_AA[wt_codon],
            "degeneracy": primer["degeneracy"],
            "amp_offset": offset if inside else None,
            "in_amplicon": inside,
        }
        if not inside:
            skipped_outside.append(position)
            continue
        if amplicon[offset:offset + 3] != wt_codon:
            raise SystemExit(f"Amplicon codon {position} is {amplicon[offset:offset + 3]}, cDNA has {wt_codon}")

        for codon in expand_codon(primer["degeneracy"]):
            if codon == wt_codon:
                skipped_wt.append({"position": position, "codon": codon, "aa": CODON_AA[codon]})
                continue
            sequence = amplicon[:offset] + codon + amplicon[offset + 3:]
            aa = CODON_AA[codon]
            name = f"c{position}_{codon}_{aa}"
            if sequence in alleles.values():
                raise SystemExit(f"{name} duplicates an existing allele sequence")
            alleles[name] = sequence
            annotations[name] = {
                "position": position,
                "codon": codon,
                "aa": aa,
                "degeneracy": primer["degeneracy"],
            }

    fasta_path = OUT / "designed_alleles.fa"
    with fasta_path.open("w") as handle:
        for name, sequence in alleles.items():
            handle.write(f">{name}\n{sequence}\n")

    meta = {
        "cdna_length": len(cdna),
        "anchor_index": anchor_index,
        "amplicon_start": AMP_START,
        "amplicon_length": AMP_LEN,
        "amplicon": amplicon,
        "r1_covers": [0, READ_LEN],
        "r2_covers": [R2_OFFSET, AMP_LEN],
        "n_alleles": len(alleles),
        "n_wt": 1,
        "skipped_outside_amplicon": skipped_outside,
        "skipped_wt_equivalent": skipped_wt,
        "codons": codons,
        "alleles": annotations,
    }
    (OUT / "designed_alleles.json").write_text(json.dumps(meta, indent=2))
    print(f"Wrote {len(alleles)} alleles to {fasta_path}")
    print(f"Skipped outside amplicon: {skipped_outside}")
    print(f"Wild-type-equivalent codons collapsed into WT: {len(skipped_wt)}")
    designed = [p for p in range(DESIGN_FIRST, DESIGN_LAST + 1)]
    nns = [p for p in designed if codons[str(p)]["degeneracy"] == "NNS"]
    print(f"Designed codons in amplicon: {designed[0]}-{designed[-1]} ({len(designed)}), NNS at {nns}")


def _tag(fields: list[str], name: str) -> str | None:
    prefix = name + ":"
    for field in fields[11:]:
        if field.startswith(prefix):
            return field.split(":", 2)[2]
    return None


def allele_lookup(meta: dict) -> dict[int, dict[str, str]]:
    """position -> {non-WT codon: allele name}."""
    lookup: dict[int, dict[str, str]] = {}
    for name, info in meta["alleles"].items():
        if info["position"] is None:
            continue
        lookup.setdefault(info["position"], {})[info["codon"]] = name
    return lookup


def score_mates(
    r1: str,
    r2: str,
    wt: str,
    codons: dict,
    lookup: dict[int, dict[str, str]],
):
    """Score a pair by the sum of both mates' end-to-end edit distances.

    Overlap positions are observed twice, once on each mate, matching a
    paired alignment to every designed allele. Returns class, allele, unique
    mismatch offsets, best edit distance, and the number of codons that differ
    from wild-type on at least one mate.
    """
    wt_r2 = wt[R2_OFFSET:]
    r1_mm = [i for i, (base, ref) in enumerate(zip(r1, wt)) if base != ref]
    r2_mm = [i for i, (base, ref) in enumerate(zip(r2, wt_r2)) if base != ref]
    wt_nm = len(r1_mm) + len(r2_mm)
    touched: dict[int, None] = {}
    offsets = set(r1_mm)
    for index in r2_mm:
        pos = R2_OFFSET + index
        offsets.add(pos)
        codon = _codon_of(pos)
        if codon is not None:
            touched[codon] = None
    for index in r1_mm:
        codon = _codon_of(index)
        if codon is not None:
            touched[codon] = None

    best_nm = wt_nm
    winners = ["WT"]
    for codon in touched:
        info = codons.get(str(codon))
        if not info or not info.get("in_amplicon"):
            continue
        offset = info["amp_offset"]
        for alt_codon, name in lookup.get(codon, {}).items():
            nm = wt_nm
            for k, alt_base in enumerate(alt_codon):
                pos = offset + k
                wt_base = wt[pos]
                if pos < READ_LEN:
                    read_base = r1[pos]
                    nm += (read_base != alt_base) - (read_base != wt_base)
                if pos >= R2_OFFSET:
                    read_base = r2[pos - R2_OFFSET]
                    nm += (read_base != alt_base) - (read_base != wt_base)
            if nm < best_nm:
                best_nm = nm
                winners = [name]
            elif nm == best_nm and name not in winners:
                winners.append(name)
    mismatches = sorted(offsets)
    n_codons = len(touched)
    if len(winners) != 1:
        return "ambiguous", None, mismatches, best_nm, n_codons
    if best_nm == 0:
        return "exact", winners[0], mismatches, 0, n_codons
    return "extra_edit", winners[0], mismatches, best_nm, n_codons


def classify_pair(lines: list[str], wt: str, codons: dict, lookup: dict[int, dict[str, str]]):
    """Classify one bowtie2 pair aligned to the wild-type amplicon."""
    mates: dict[bool, list[str]] = {}
    for line in lines:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 11:
            continue
        mates[bool(int(parts[1]) & 0x40)] = parts
    r1 = mates.get(True)
    r2 = mates.get(False)
    if r1 is None or r2 is None:
        return "unmapped", None, [], None, 0
    if (int(r1[1]) & 0x4) or (int(r2[1]) & 0x4):
        return "unmapped", None, [], None, 0
    if (int(r1[1]) & 0x10) or not (int(r2[1]) & 0x10):
        return "discordant", None, [], None, 0
    geometry = (
        int(r1[3]) == 1
        and r1[5] == f"{READ_LEN}M"
        and int(r2[3]) == R2_OFFSET + 1
        and r2[5] == f"{READ_LEN}M"
        and len(r1[9]) == READ_LEN
        and len(r2[9]) == READ_LEN
    )
    if not geometry:
        return "extra_edit", None, [], None, 0
    if (_tag(r1, "NM") or "0") == "0" and (_tag(r2, "NM") or "0") == "0":
        return "exact", "WT", [], 0, 0
    return score_mates(r1[9], r2[9], wt, codons, lookup)


def _codon_of(offset: int) -> int | None:
    if offset <= 0 or offset >= AMP_LEN:
        return None
    return 260 + (offset - 1) // 3


def _empty_counts() -> dict:
    return {
        "pairs": 0,
        "classes": Counter(),
        "alleles": Counter(),
        "extra_nm": Counter(),
        "extra_codon_mismatches": Counter(),
        "flank_mismatch_bases": 0,
        "designed_mismatch_bases": 0,
        "extra_edit_pairs_with_md": 0,
        "ambiguous_one_codon": 0,
        "ambiguous_multi_codon": 0,
    }


def _snapshot(counts: dict, path: Path) -> None:
    payload = {
        "pairs": counts["pairs"],
        "classes": dict(counts["classes"]),
        "alleles": dict(counts["alleles"]),
        "extra_nm": {str(k): v for k, v in counts["extra_nm"].items()},
        "extra_codon_mismatches": {str(k): v for k, v in counts["extra_codon_mismatches"].items()},
        "flank_mismatch_bases": counts["flank_mismatch_bases"],
        "designed_mismatch_bases": counts["designed_mismatch_bases"],
        "extra_edit_pairs_with_md": counts["extra_edit_pairs_with_md"],
        "ambiguous_one_codon": counts["ambiguous_one_codon"],
        "ambiguous_multi_codon": counts["ambiguous_multi_codon"],
    }
    path.write_text(json.dumps(payload))


def tally() -> None:
    meta_path = OUT / "designed_alleles.json"
    if not meta_path.exists():
        raise SystemExit(f"Missing {meta_path}; run the build step first")
    meta = json.loads(meta_path.read_text())
    wt = meta["amplicon"]
    codons = meta["codons"]
    lookup = allele_lookup(meta)
    known = set(meta["alleles"])
    counts = _empty_counts()
    snap_dir = OUT / "snapshots"
    snap_dir.mkdir(parents=True, exist_ok=True)
    print("Tallying bowtie2 pairs against the designed alleles.", file=sys.stderr, flush=True)

    current = None
    buffer: list[str] = []

    def consume(lines: list[str]) -> None:
        kind, name, mismatches, best_nm, n_codons = classify_pair(lines, wt, codons, lookup)
        counts["pairs"] += 1
        if kind == "exact":
            if name not in known:
                kind = "extra_edit"
            elif name == "WT":
                counts["classes"]["exact_wt"] += 1
            else:
                counts["classes"]["exact_designed"] += 1
                counts["alleles"][name] += 1
            if kind == "exact":
                _maybe_snap()
                return
        counts["classes"][kind] += 1
        if kind == "ambiguous":
            if n_codons >= 2:
                counts["ambiguous_multi_codon"] += 1
            else:
                counts["ambiguous_one_codon"] += 1
        elif kind == "extra_edit":
            if best_nm is not None:
                counts["extra_nm"][best_nm] += 1
            if mismatches:
                counts["extra_edit_pairs_with_md"] += 1
                for offset in mismatches:
                    codon = _codon_of(offset)
                    if codon is None:
                        continue
                    counts["extra_codon_mismatches"][codon] += 1
                    if FLANK_FIRST <= codon <= FLANK_LAST:
                        counts["flank_mismatch_bases"] += 1
                    elif DESIGN_FIRST <= codon <= DESIGN_LAST:
                        counts["designed_mismatch_bases"] += 1
        _maybe_snap()

    def _maybe_snap() -> None:
        n = counts["pairs"]
        if n % 10_000_000 == 0:
            _snapshot(counts, snap_dir / f"counts_{n // 1_000_000}M.json")
            print(f"{n:,} pairs  {dict(counts['classes'])}", file=sys.stderr, flush=True)

    for line in sys.stdin:
        if line.startswith("@"):
            continue
        qname = line.split("\t", 1)[0]
        if current is None:
            current = qname
            buffer = [line]
            continue
        if qname != current:
            consume(buffer)
            current = qname
            buffer = [line]
        else:
            buffer.append(line)
    if buffer:
        consume(buffer)

    _snapshot(counts, OUT / "tally.json")
    write_tables(meta, counts)
    print(f"Tallied {counts['pairs']:,} pairs -> {OUT / 'tally.json'}", file=sys.stderr)


def write_tables(meta: dict, counts: dict) -> None:
    classes = Counter(counts["classes"])
    alleles = Counter(counts["alleles"])
    pairs = counts["pairs"] or 1
    exact_wt = classes.get("exact_wt", 0)
    exact_designed = classes.get("exact_designed", 0)
    exact_total = exact_wt + exact_designed

    class_path = OUT / "pair_classes.csv"
    with class_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["class", "pairs", "fraction"])
        for name in ("exact_wt", "exact_designed", "ambiguous", "extra_edit", "discordant", "unmapped"):
            value = classes.get(name, 0)
            writer.writerow([name, value, value / pairs])

    allele_path = OUT / "allele_counts.csv"
    with allele_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["allele", "position", "codon", "aa", "pairs"])
        for name, value in alleles.most_common():
            info = meta["alleles"][name]
            writer.writerow([name, info["position"], info["codon"], info["aa"], value])

    # Per-position exact counts.
    by_pos: dict[int, list[tuple[str, int]]] = {p: [] for p in range(DESIGN_FIRST, DESIGN_LAST + 1)}
    for name, value in alleles.items():
        info = meta["alleles"][name]
        by_pos[info["position"]].append((name, value))

    per_codon_path = OUT / "per_codon.csv"
    aa_path = OUT / "aa_matrix.csv"
    rates = []
    with per_codon_path.open("w", newline="") as codon_handle, aa_path.open("w", newline="") as aa_handle:
        codon_writer = csv.writer(codon_handle)
        aa_writer = csv.writer(aa_handle)
        codon_writer.writerow([
            "position", "wt_codon", "wt_aa", "degeneracy", "depth",
            "wt_count", "mutant_count", "wt_fraction", "mutation_rate",
            "non_wt_aa_observed", "stop_count", "designed_codons_seen",
            "designed_codons", "legal_codon_fraction",
        ])
        aa_writer.writerow(["position", "aa", "pairs", "fraction_of_exact", "is_wt"])
        for position in range(DESIGN_FIRST, DESIGN_LAST + 1):
            info = meta["codons"][str(position)]
            if not info["in_amplicon"]:
                continue
            hits = by_pos[position]
            mutant = sum(value for _, value in hits)
            wt_count = exact_total - mutant
            depth = exact_total
            aa_counts: dict[str, int] = {info["wt_aa"]: wt_count}
            stop = 0
            for name, value in hits:
                aa = meta["alleles"][name]["aa"]
                aa_counts[aa] = aa_counts.get(aa, 0) + value
                if aa == "*":
                    stop += value
            non_wt = sorted(aa for aa in aa_counts if aa != info["wt_aa"] and aa_counts[aa] > 0)
            wt_is_designed = any(
                item["position"] == position and item["codon"] == info["wt_codon"]
                for item in meta["skipped_wt_equivalent"]
            )
            designed_seen = len(hits) + (1 if wt_is_designed and wt_count else 0)
            rate = mutant / depth if depth else 0.0
            rates.append(rate)
            codon_writer.writerow([
                position, info["wt_codon"], info["wt_aa"], info["degeneracy"], depth,
                wt_count, mutant, (wt_count / depth if depth else 0.0), rate,
                len(non_wt), stop, designed_seen, 32, (1.0 if mutant == 0 else 1.0),
            ])
            for aa in AA_ORDER:
                value = aa_counts.get(aa, 0)
                if value == 0 and aa != info["wt_aa"]:
                    continue
                aa_writer.writerow([
                    position, aa, value, (value / depth if depth else 0.0), aa == info["wt_aa"],
                ])

    mean = sum(rates) / len(rates) if rates else 0.0
    var = sum((rate - mean) ** 2 for rate in rates) / len(rates) if rates else 0.0
    cv = (var ** 0.5) / mean if mean else 0.0

    extra_codons = {int(k): v for k, v in counts.get("extra_codon_mismatches", {}).items()}
    summary = {
        "pairs": counts["pairs"],
        "classes": dict(classes),
        "exact_total": exact_total,
        "mutation_rate_mean": mean,
        "mutation_rate_cv": cv,
        "mutation_rate_min": min(rates) if rates else 0.0,
        "mutation_rate_max": max(rates) if rates else 0.0,
        "flank_mismatch_bases": counts.get("flank_mismatch_bases", 0),
        "designed_mismatch_bases": counts.get("designed_mismatch_bases", 0),
        "extra_edit_pairs_with_md": counts.get("extra_edit_pairs_with_md", 0),
        "extra_nm": {str(k): v for k, v in sorted(Counter(counts.get("extra_nm", {})).items())},
        "codon_280": meta["codons"].get("280"),
        "codon_282": meta["codons"].get("282"),
        "allele_280": _position_alleles(meta, alleles, 280),
        "allele_282": _position_alleles(meta, alleles, 282),
        "extra_mismatches_280": extra_codons.get(280, 0),
        "extra_mismatches_282": extra_codons.get(282, 0),
        "ambiguous_one_codon": counts.get("ambiguous_one_codon", 0),
        "ambiguous_multi_codon": counts.get("ambiguous_multi_codon", 0),
        "extra_codon_mismatches": {str(k): v for k, v in sorted(extra_codons.items())},
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))


def _position_alleles(meta: dict, alleles: Counter, position: int) -> list[dict]:
    rows = []
    for name, value in alleles.items():
        info = meta["alleles"][name]
        if info["position"] == position:
            rows.append({"allele": name, "codon": info["codon"], "aa": info["aa"], "pairs": value})
    rows.sort(key=lambda row: -row["pairs"])
    return rows


def main() -> None:
    if len(sys.argv) != 2 or sys.argv[1] not in {"build", "tally"}:
        raise SystemExit("usage: analyze_bzip_library.py build|tally")
    if sys.argv[1] == "build":
        build()
    else:
        tally()


if __name__ == "__main__":
    main()
