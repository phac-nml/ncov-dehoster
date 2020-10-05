process PLACEHOLDER {
    publishDir "${params.outdir}/test_dir", pattern: "*", mode: "copy"

    input:
    tuple(sampleName, path(forward), path(reverse))

    output:
    path("*.fastq")
    path("*.txt")

    script:
    """
    touch ${sampleName}.txt
    mv ${forward} forward_${forward}
    mv ${reverse} forward_${reverse}
    """
}

process captureViralReads {

    publishDir "${params.outdir}/virusSAM", pattern: "*.virus.sorted.sam", mode: "copy"

    label 'largecpu'

    input:
    tuple(sampleName, path(forward), path(reverse), path(covid_reference))

    output:
    tuple(sampleName, path("*.virus.sorted.sam"))

    script:
    """
    bwa index ${covid_reference}
    bwa mem -t 16 ${covid_reference} ${forward} ${reverse} | samtools sort --threads 6 -T "temp" -O SAM -o ${sampleName}.virus.sorted.sam
    """
}

process captureHumanReads {

    publishDir "${params.outdir}/humanSAM", pattern: "*.human.sorted.bam", mode: "copy"

    label 'bwa_human'

    input:
    tuple(sampleName, path(forward), path(reverse), path(human_reference))


    output:
    path("*.human.sorted.bam")

    script:
    """
    bwa index ${human_reference}
    bwa mem -t 16 ${human_reference} ${forward} ${reverse} | samtools sort --threads 10 -T "temp" -O SAM -o ${sampleName}.human.sorted.sam
    """
}