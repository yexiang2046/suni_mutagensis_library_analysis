# SUNi library count

Counts a saturation-mutagenesis library the way the SUNi paper does. Read pairs are merged, merges with expected error of 0.5 or more are dropped, a 3' adapter is trimmed when present, and a read is kept only when it is a perfect match to one designed `NNK` or `NNS` allele.

Method: Mighell, Toledano, and Lehner, *PLOS ONE* 2023. [10.1371/journal.pone.0288158](https://doi.org/10.1371/journal.pone.0288158)

## Run

Needs `nextflow`, `vsearch`, `cutadapt`, and `python3` on `PATH`. The count scripts use only the Python standard library. Read pairs are inputs and are not part of this repository. An example cDNA sequence and primer table are in `testdata/`. The cDNA is a FASTA or plain DNA file.

```bash
nextflow run main.nf \
    --read1 R1.fastq.gz \
    --read2 R2.fastq.gz \
    --cdna template.fa \
    --primers primers.csv
```

`build_alleles` writes `designed_alleles.fa` and `designed_alleles.json`. `count_library` streams the merge, the adapter trim, and the perfect-match count, so the merged FASTQ is not stored. The workflow publishes five files to `results/`:

- `designed_alleles.fa`
- `designed_alleles.json`
- `paper_style_summary.json`
- `vsearch.log`
- `cutadapt.log`

| Parameter | Role |
| --- | --- |
| `anchor` | First 30 bp of the ORF, starting at the start codon |
| `amp_start`, `amp_len` | Amplicon coordinates on the cDNA |
| `design_first`, `design_last` | Inclusive codon range in the primer table |
| `maxee` | Maximum expected error of a merged read |
| `adapter`, `adapter_overlap` | 3' adapter and the minimum overlap to trim |
| `threads` | Threads for the merge |

A designed codon that falls outside the amplicon is left out of the reference. Extending `amp_len` is how a later codon enters the count. `nextflow.config` supplies one example of these values.

## Scripts

`analysis_library.py build` constructs the allele reference from the cDNA sequence and primer CSV. `count.py` reads merged FASTQ on stdin and writes the summary. Those two commands are the processes in `main.nf`.

```bash
python3 analysis_library.py build \
    --cdna template.fa \
    --primers primers.csv \
    --outdir alleles

vsearch --fastq_mergepairs R1.fastq.gz --reverse R2.fastq.gz \
    --fastq_maxee 0.5 --threads 8 --fastqout - \
| cutadapt -a AGATCGGAAGAGC --overlap 10 -o - - \
| python3 count.py --meta alleles/designed_alleles.json summary.json
```

The summary reports retained reads, wild-type, programmed, and other fractions, per-codon frequencies, LogDiff, and screening efficiency. LogDiff is log10 of the 90th/10th percentile of per-codon programmed frequency. Screening efficiency is the programmed fraction divided by 10^LogDiff.
