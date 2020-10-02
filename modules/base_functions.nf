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

    publishDir "${params.outdir}/output_fastq", pattern: "${sampleName}.cleaned.fastq", mode: "copy"

    label 'smallcpu'

    tag { sampleName }

    input:
    tuple sampleName, file(unmapped_bam)

    output:
    file("${sampleName}.cleaned.fastq")

    script:
    """
    samtools fastq --threads 4 ${unmapped_bam} > ${sampleName}.cleaned.fastq
    """
}
