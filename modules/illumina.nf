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

    label 'largecpu'

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

    publishDir "${params.outdir}/virusSAM", pattern: "${sampleName}.*", mode: "copy"

    label 'largecpu'

    input:
    tuple(sampleName, path(forward), path(reverse), path(covid_reference))
    path(indexed_reference)

    output:
    tuple sampleName, path("${sampleName}.virus.sorted.sam"), emit: sam
    path("${sampleName}.flagstats.txt")

    script:
    """
    bwa mem -t 10 ${covid_reference} ${forward} ${reverse} | samtools sort --threads 6 -T "temp" -O SAM -o ${sampleName}.sorted.sam
    samtools view --threads 4 -h -f 0x0002 ${sampleName}.sorted.sam > ${sampleName}.virus.sorted.sam
    samtools flagstat ${sampleName}.virus.sorted.sam > ${sampleName}.flagstats.txt
    """
}

process captureHumanReads {

    publishDir "${params.outdir}/humanSAM", pattern: "${sampleName}.*", mode: "copy"

    label 'largecpu'

    input:
    tuple(sampleName, path(forward), path(reverse), path(human_reference))
    path(indexed_reference)

    output:
    tuple sampleName, path("${sampleName}.human.sorted.sam"), emit: sam
    path("${sampleName}.flagstats.txt")

    script:
    """
    bwa mem -t 10 ${human_reference} ${forward} ${reverse} | samtools sort --threads 6 -T "temp" -O SAM -o ${sampleName}.sorted.sam
    samtools view --threads 4 -h -F 4 ${sampleName}.sorted.sam > ${sampleName}.human.sorted.sam
    samtools flagstat ${sampleName}.human.sorted.sam > ${sampleName}.flagstats.txt
    """
}

process dehostSamFiles {

    publishDir "${params.outdir}/dehosted_sams", pattern: "${sampleName}.*", mode: "copy"

    label 'smallcpu'

    input:
    tuple(sampleName, path(virus_sam))
    tuple(sampleName, path(human_sam))

    output:
    tuple sampleName, path("${sampleName}.dehosted.sam")

    script:
    """
    dehost.py --keep_sam ${virus_sam} --remove_sam ${human_sam} -q ${params.keep_min_map_quality} -Q ${params.remove_min_map_quality} -o ${sampleName}.dehosted.sam
    """
}

process generateDehostedReads {

    publishDir "${params.outdir}/dehosted_paired_fastqs", pattern: "${sampleName}_R*", mode: "copy"

    input:
    tuple(sampleName, path(dehosted_sam))

    output:
    path("${sampleName}_R*")

    script:
    """
    samtools fastq -1 ${sampleName}_R1.dehosted.fastq -2 ${sampleName}_R2.dehosted.fastq ${dehosted_sam}
    """
}
