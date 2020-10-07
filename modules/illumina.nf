process indexHumanReference {

    label 'bwa_human_index'

    input:
    path(human_ref)
    path(bwa_index_folder)

    output:
    file("*")

    script:
    if ( params.human_bwa_index )

        """
        ln -s $bwa_index_folder/* .
        """

    else

        """
        bwa index -a bwtsw $human_ref
        """
}

process indexViralReference {

    label 'mediumcpu'

    input:
    path(viral_ref)

    output:
    file("*")

    script:
    """
    bwa index $viral_ref
    """
}

process captureViralReads {

    publishDir "${params.outdir}/virusMAPs", pattern: "${sampleName}.*", mode: "copy"

    label 'largecpu'

    input:
    tuple(sampleName, path(forward), path(reverse), path(covid_reference))
    path(indexed_reference)

    output:
    tuple sampleName, path("${sampleName}.virus.sorted.bam"), emit: bam
    path("${sampleName}.flagstats.txt")

    script:
    """
    bwa mem -t 10 ${covid_reference} ${forward} ${reverse} | samtools sort --threads 6 -T "temp" -O BAM -o ${sampleName}.sorted.bam
    samtools view --threads 4 -h -f 0x0002 ${sampleName}.sorted.bam > ${sampleName}.virus.sorted.bam
    samtools flagstat ${sampleName}.virus.sorted.bam > ${sampleName}.flagstats.txt
    """
}

process captureHumanReads {

    publishDir "${params.outdir}/humanMAPs", pattern: "${sampleName}.*", mode: "copy"

    label 'largecpu'

    input:
    tuple(sampleName, path(forward), path(reverse), path(human_reference))
    path(indexed_reference)

    output:
    tuple sampleName, path("${sampleName}.human.sorted.bam"), emit: bam
    path("${sampleName}.flagstats.txt")

    script:
    """
    bwa mem -t 10 ${human_reference} ${forward} ${reverse} | samtools sort --threads 6 -T "temp" -O BAM -o ${sampleName}.sorted.bam
    samtools view --threads 4 -h -F 4 ${sampleName}.sorted.bam > ${sampleName}.human.sorted.bam
    samtools flagstat ${sampleName}.human.sorted.bam > ${sampleName}.flagstats.txt
    """
}

process dehostBamFiles {

    publishDir "${params.outdir}/dehostedBAMs", pattern: "${sampleName}.dehosted.bam", mode: "copy"

    label 'mediumcpu'

    input:
    tuple(sampleName, path(virus_bam), path(human_bam))

    output:
    tuple sampleName, path("${sampleName}.dehosted.bam"), emit: bam
    path "${sampleName}*.csv", emit: csv

    script:
    """
    dehost.py --keep ${virus_bam} --remove ${human_bam} -q ${params.keep_min_map_quality} -Q ${params.remove_min_map_quality} -o ${sampleName}.dehosted.bam
    """
}

process generateDehostedReads {

    publishDir "${params.outdir}/dehosted_paired_fastqs", pattern: "${sampleName}-dehosted_R*", mode: "copy"

    label 'mediumcpu'

    input:
    tuple(sampleName, path(dehosted_bam))

    output:
    path("${sampleName}-dehosted_R*")

    script:
    """
    samtools fastq -1 ${sampleName}-dehosted_R1.fastq -2 ${sampleName}-dehosted_R2.fastq ${dehosted_bam}
    """
}

process combineCSVs {

    publishDir "${params.outdir}", pattern: "removal_summary.csv", mode: "copy"

    label 'smallcpu'

    input:
    path(csvs)

    output:
    path("removal_summary.csv")

    script:
    """
    csvtk concat *.virus_stats.csv > removal_summary.csv
    """
}
