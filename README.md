# JunB bZIP SUNi library analysis

Counts a saturation-mutagenesis library of the JunB bZIP the way the SUNi paper does: merge read pairs, drop merges with expected error of 0.5 or more, trim a 3' TruSeq adapter, then keep a read only when it is a perfect match to one designed `NNK` or `NNS` allele.

Method: Mighell, Toledano, and Lehner, *PLOS ONE* 2023. [10.1371/journal.pone.0288158](https://doi.org/10.1371/journal.pone.0288158)

Sample `15485-RB-0001` is a 211 bp amplicon (cDNA 1062–1272) sequenced as 151 bp paired reads. Designed mutagenesis covers codons 267–330. Codon 330 is the first base after the amplicon, so it is omitted from the reference. A longer amplicon is how it would be included later.

## Result

Of 176,134,252 pairs, 164,400,394 (93.3%) were retained.

| Class | Reads | Fraction of retained |
| --- | ---: | ---: |
| Wild type | 142,427,245 | 86.63% |
| Programmed `NNK`/`NNS` | 9,491,738 | 5.77% |
| Other | 12,481,411 | 7.59% |

LogDiff, log10 of the 90th/10th percentile of per-codon programmed frequency, is 1.41. Screening efficiency (programmed fraction / 10^LogDiff) is 0.0023. The full table is `bzip_mutation_profile/paper_style_summary.json`.

## Run

Needs `nextflow`, `vsearch`, `cutadapt`, and `python3` on `PATH`. The count scripts use only the Python standard library.

The read pairs, `JunB cDNA.docx`, and the bZIP primer CSV are not in this repository. Place them at the paths in `nextflow.config`, or pass them on the command line:

```bash
nextflow run main.nf \
    --read1 R1.fastq.gz \
    --read2 R2.fastq.gz \
    --cdna "JunB cDNA.docx" \
    --primers "JunB SUNI Primers/primers_output_JunB_bZIP_codons267-330.csv"
```

`build_alleles` writes `designed_alleles.fa` and `designed_alleles.json`. `count_library` streams the merge, the adapter trim, and the perfect-match count, so the merged FASTQ is not stored. Published files go to `results/`: the summary plus `vsearch.log` and `cutadapt.log`.

Defaults in `nextflow.config` are this library: amplicon start 1062, length 211, designed codons 267–330, expected error 0.5, adapter `AGATCGGAAGAGC`, overlap 10, 8 threads.

## Scripts

`analyze_bzip_library.py build` constructs the allele reference from the cDNA and primer CSV. `tally` is a separate bowtie2 end-to-end count against the wild-type amplicon. It is not a Nextflow process. The paper-style numbers above come from `count_paper_style.py`.

```bash
python3 analyze_bzip_library.py build --outdir bzip_mutation_profile
vsearch --fastq_mergepairs R1.fastq.gz --reverse R2.fastq.gz \
    --fastq_maxee 0.5 --threads 8 --fastqout - \
| cutadapt -a AGATCGGAAGAGC --overlap 10 -o - - \
| python3 count_paper_style.py \
    --meta bzip_mutation_profile/designed_alleles.json \
    paper_style_summary.json
```
