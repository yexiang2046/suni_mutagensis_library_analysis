#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Paper-style SUNi library count.
 * vsearch merges pairs and drops expected error >= maxee, cutadapt trims a
 * 3' adapter when present, and count_paper_style.py keeps only a perfect
 * match to one designed NNK/NNS allele. The merged FASTQ is streamed.
 */

process build_alleles {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path cdna
    path primers

    output:
    path "designed_alleles.fa", emit: fasta
    path "designed_alleles.json", emit: json

    script:
    """
    python3 ${projectDir}/analyze_bzip_library.py build \
        --cdna ${cdna} \
        --primers ${primers} \
        --outdir . \
        --anchor ${params.anchor} \
        --amp-start ${params.amp_start} \
        --amp-len ${params.amp_len} \
        --design-first ${params.design_first} \
        --design-last ${params.design_last}
    """
}

process count_library {
    cpus params.threads
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path read1
    path read2
    path alleles_json

    output:
    path "paper_style_summary.json"
    path "vsearch.log"
    path "cutadapt.log"

    script:
    """
    set -o pipefail
    vsearch --fastq_mergepairs ${read1} --reverse ${read2} \\
        --fastq_maxee ${params.maxee} --threads ${task.cpus} --fastqout - \\
        2> vsearch.log \\
    | cutadapt -a ${params.adapter} --overlap ${params.adapter_overlap} -o - - \\
        2> cutadapt.log \\
    | python3 ${projectDir}/count_paper_style.py --meta ${alleles_json} paper_style_summary.json
    """
}

workflow {
    alleles = build_alleles(
        file(params.cdna, checkIfExists: true),
        file(params.primers, checkIfExists: true),
    )
    count_library(
        file(params.read1, checkIfExists: true),
        file(params.read2, checkIfExists: true),
        alleles.json,
    )
}
