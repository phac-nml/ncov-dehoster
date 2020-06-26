process copyReference {

    label 'smallcpu'

    input:
    file(reference)

    output:
    file('reference.fasta')
    

    script:
    """
    mv ${reference} ./reference.fasta
    """
}

process minimap2 {

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.sorted.bam", mode: "copy"

    label 'minimap'

    tag { sampleName }

    input:
    tuple path(fastq), file(reference)

    output:
    tuple sampleName, file("${sampleName}.sorted.bam")
    

    script:
    sampleName = fastq.getBaseName().replaceAll(~/\.fastq.*$/, '')

    """
    minimap2 -x map-ont -t 160 ${reference} ${fastq} -a | samtools sort --threads 10 -T "temp" -O BAM -o ${sampleName}.sorted.bam
    """
}

process removeMappedReads {

    label 'smallcpu'

    tag { sampleName }

    input:
    tuple sampleName, file(sorted_bam)

    output:
    tuple sampleName, file("${sampleName}.unmapped.sorted.bam")
    

    script:
    """
    samtools view -f4 -b --threads 4 -o ${sampleName}.unmapped.sorted.bam ${sorted_bam}
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