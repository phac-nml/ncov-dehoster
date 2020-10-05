process copyReference {

    label 'smallcpu'

    input:
    file(human_ref)

    output:
    file('reference.fasta')
    

    script:
    """
    mv ${human_ref} ./reference.fasta
    """
}

process extractFastq {

    publishDir "${params.outdir}/output_fastq", pattern: "${sampleName}.mapped.fastq", mode: "copy"

    label 'smallcpu'

    tag { sampleName }

    input:
    tuple sampleName, path(mapped_bam)

    output:
    file("${sampleName}.mapped.fastq")

    script:
    """
    samtools fastq --threads 4 ${mapped_bam} > ${sampleName}.mapped.fastq
    """
}
