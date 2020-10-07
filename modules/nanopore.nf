process minimap2 {

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.sorted.bam", mode: "copy"

    label 'minimap'

    tag { sampleName }

    input:
    tuple path(fastq), file(human_ref), file(cov2019_ref)

    output:
    tuple sampleName, file("${sampleName}.sorted.bam")
    

    script:
    sampleName = fastq.getBaseName().replaceAll(~/\.fastq.*$/, '')

    """
    minimap2 -x map-ont -t 160 ${human_ref} ${cov2019_ref} ${fastq} -a | samtools sort --threads 10 -T "temp" -O BAM -o ${sampleName}.sorted.bam
    """
}

process samtoolsFlagstat {

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.flagstats.txt", mode: "copy"

    label 'smallcpu'

    tag { sampleName }

    input:
    tuple sampleName, file(sorted_bam)

    output:

    file("${sampleName}.flagstats.txt") 

    script:
    """
    samtools flagstat ${sorted_bam} > ${sampleName}.flagstats.txt
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
