process generateCompositeReference {
    label 'smallCPU'

    input:
    path(human_ref)
    path(viral_ref)

    output:
    file("composite_reference.fa")

    script:
    """
    cat $human_ref $viral_ref > composite_reference.fa
    """
}

process grabCompositeIndex {
    label 'smallCPU'

    input:
    path(index_folder)

    output:
    file("*.fa*")

    script:
    """
    ln -sf $index_folder/*.fa* ./
    """
}

process indexCompositeReference {
    publishDir "${params.outdir}/humanBWAIndex", pattern: "*.fa*", mode: "symlink"

    label 'indexResources'

    input:
    path(composite_ref)

    output:
    path("*.fa*"), emit: index
    path("*.process.yml"), emit: versions

    script:
    """
    bwa index -a bwtsw $composite_ref

    # Versions #
    cat <<-END_VERSIONS > indexbwa.process.yml
        "${task.process}":
            bwa: \$(echo \$(bwa 2> >(grep -o 'Version' | cut -d ' ' -f 2)))
    END_VERSIONS
    """
}

process compositeMappingBWA {
    publishDir "${params.outdir}/compositeMAPs", pattern: "${sampleName}.*", mode: "copy"

    input:
    tuple( val(sampleName), path(forward), path(reverse), path(composite_reference))
    path(indexed_reference)

    output:
    tuple val(sampleName), path("${sampleName}.sorted.bam"), emit: bam
    path("${sampleName}.flagstats.txt"), emit: flagstats
    path("*.process.yml"), emit: versions

    script:
    processThreads = task.cpus * 2
    """
    bwa mem -t ${processThreads} ${composite_reference} ${forward} ${reverse} | samtools sort --threads ${task.cpus} -T "temp" -O BAM -o ${sampleName}.sorted.bam
    samtools flagstat ${sampleName}.sorted.bam > ${sampleName}.flagstats.txt

    # Versions #
    cat <<-END_VERSIONS > mapbwa.process.yml
        "${task.process}":
            bwa: \$(echo \$(bwa 2> >(grep -o 'Version' | cut -d ' ' -f 2)))
            samtools: \$(echo \$(samtools --version | head -n 1 | grep samtools | sed 's/samtools //'))
    END_VERSIONS
    """
}

process dehostBamFiles {
    publishDir "${params.outdir}/dehostedBAMs", pattern: "${sampleName}.dehosted.bam", mode: "copy"

    label 'smallCPU'
    tag { sampleName }

    input:
    tuple( val(sampleName), path(composite_bam))

    output:
    tuple val(sampleName), path("${sampleName}.dehosted.bam"), emit: bam
    path "${sampleName}*.csv", emit: csv
    path("*.process.yml"), emit: versions

    script:

    def rev = workflow.commitId ?: workflow.revision ?: workflow.scriptId

    """
    samtools index ${composite_bam}
    dehost_illumina.py --file ${composite_bam} \
    --keep_id ${params.keep_ref_id} \
    -q ${params.keep_min_map_quality} \
    -Q ${params.remove_min_map_quality} \
    -o ${sampleName}.dehosted.bam \
    -R ${rev} 

    # Versions #
    cat <<-END_VERSIONS > dehostbam.process.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version | head -n 1 | grep samtools | sed 's/samtools //'))
            dehost_illumina.py: 0.1.0
    END_VERSIONS
    """
}

process generateDehostedReads {
    publishDir "${params.outdir}/dehosted_paired_fastqs", pattern: "${sampleName}_dehosted_R*", mode: "copy"

    label 'mediumMem'
    tag { sampleName }

    input:
    tuple( val(sampleName), path(dehosted_bam))

    output:
    path("${sampleName}_dehosted_R*"), emit: dehosted_fastq
    path("*.process.yml"), emit: versions

    script:
    """
    samtools fastq -1 ${sampleName}_dehosted_R1.fastq.gz -2 ${sampleName}_dehosted_R2.fastq.gz ${dehosted_bam}

    # Versions #
    cat <<-END_VERSIONS > index.process.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version | head -n 1 | grep samtools | sed 's/samtools //'))
    END_VERSIONS
    """
}

process combineCSVs {
    publishDir "${params.outdir}", pattern: "removal_summary.csv", mode: "copy"

    label 'smallCPU'

    input:
    path(csvs)

    output:
    path("removal_summary.csv"), emit: summary
    path("*.process.yml"), emit: versions

    script:
    """
    csvtk concat *_stats.csv > removal_summary.csv

    # Versions #
    cat <<-END_VERSIONS > csv.process.yml
        "${task.process}":
            csvtk: \$(echo \$(csvtk version | sed 's/csvtk //'))
    END_VERSIONS
    """
}
