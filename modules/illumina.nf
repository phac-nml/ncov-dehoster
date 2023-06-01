process generateCompositeReference {
    label 'smallCPU'

    input:
    path(human_ref)
    path(viral_ref)

    output:
    path("composite_reference.fa")

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
    path("*.fa*")

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
    path("versions.yml"), emit: versions

    script:
    """
    bwa index -a bwtsw $composite_ref

    # Versions #
    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bwa: \$(echo \$(bwa 2> >(grep -o 'Version' | cut -d ' ' -f 2)))
    END_VERSIONS
    """
}

process compositeMappingBWA {
    publishDir "${params.outdir}/compositeMAPs", pattern: "${sampleName}.*", mode: "copy"

    input:
    tuple val(sampleName), path(forward), path(reverse), path(composite_reference)
    path(indexed_reference)

    output:
    tuple val(sampleName), path("${sampleName}.sorted.bam"), emit: bam
    path("${sampleName}.flagstats.txt"), emit: flagstats
    path("versions.yml"), emit: versions

    script:
    processThreads = task.cpus * 2
    """
    bwa mem -t ${processThreads} ${composite_reference} ${forward} ${reverse} \\
        | samtools sort --threads ${task.cpus} -T "temp" -O BAM -o ${sampleName}.sorted.bam
    samtools flagstat ${sampleName}.sorted.bam > ${sampleName}.flagstats.txt

    # Versions #
    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bwa: \$(echo \$(bwa 2> >(grep -o 'Version' | cut -d ' ' -f 2)))
            samtools: \$(echo \$(samtools --version | head -n 1 | grep samtools | sed 's/samtools //'))
    END_VERSIONS
    """
}

process dehostBamFiles {
    publishDir "${params.outdir}/dehostedBAMs", pattern: "${sampleName}.dehosted.sorted.bam", mode: "copy"

    label 'smallCPU'
    tag { sampleName }

    input:
    tuple val(sampleName), path(composite_bam)

    output:
    tuple val(sampleName), path("${sampleName}.dehosted.sorted.bam"), emit: bam
    path("${sampleName}*.csv"), emit: csv
    path("versions.yml"), emit: versions

    script:
    def rev = workflow.commitId ?: workflow.revision ?: workflow.scriptId
    downsample_args = []
    if ( params.downsample ) {
        downsample_args.add("--downsampled") 
        downsample_args.add("--downsampled_count $params.downsample_count")
        downsample_args.add("--downsampled_seed $params.downsample_seed")
    }
    downsample_args_final = downsample_args.join(" ")
    """
    samtools index ${composite_bam}
    dehost_illumina.py \\
        $downsample_args_final \\
        --file ${composite_bam} \\
        --keep_id ${params.keep_ref_id} \\
        -q ${params.keep_min_map_quality} \\
        -Q ${params.remove_min_map_quality} \\
        -o ${sampleName}.dehosted.bam \\
        -R ${rev}
    samtools sort ${sampleName}.dehosted.bam > ${sampleName}.dehosted.sorted.bam

    # Versions #
    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version | head -n 1 | grep samtools | sed 's/samtools //'))
            dehost_illumina.py: 0.1.0
    END_VERSIONS
    """
}

process generateDehostedReads {
    // Output if we are not fastq downsampling
    if ( !params.downsample || params.downsample_amplicons ) {
        publishDir "${params.outdir}/dehosted_paired_fastqs", pattern: "*_dehosted_R*", mode: "copy"
    }

    label 'mediumMem'
    tag { sampleName }

    input:
    tuple val(sampleName), path(dehosted_bam)

    output:
    tuple val(sampleName), path("${sampleName}_dehosted_R*"), emit: dehosted_fastq
    path("versions.yml"), emit: versions

    script:
    """
    samtools fastq -1 ${sampleName}_dehosted_R1.fastq.gz -2 ${sampleName}_dehosted_R2.fastq.gz ${dehosted_bam}

    # Versions #
    cat <<-END_VERSIONS > versions.yml
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
    path("versions.yml"), emit: versions

    script:
    """
    csvtk concat *_stats.csv > removal_summary.csv

    # Versions #
    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            csvtk: \$(echo \$(csvtk version | sed 's/csvtk //'))
    END_VERSIONS
    """
}
