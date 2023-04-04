process generateMinimap2Index {
    publishDir "${params.outdir}/compositeMinimapIndex", pattern: "composite_ref.mmi", mode: "symlink"

    label 'indexResources'

    input:
    path(human_ref)
    path(viral_ref)

    output:
    path("composite_ref.mmi"), emit: index
    path("*.process.yml"), emit: versions

    script:
    """
    cat $human_ref $viral_ref > composite_reference.fa
    minimap2 -d composite_ref.mmi composite_reference.fa

    # Versions #
    cat <<-END_VERSIONS > index.process.yml
        "${task.process}":
            minimap2: \$(echo \$(minimap2 --version))
    END_VERSIONS
    """
}

process fastqSizeSelection_MM2 {
    label 'mediumMem'
    tag { sampleName }

    input:
    path(fastq)

    output:
    tuple val(sampleName), path("${sampleName}_size_selected*.fastq"), emit: fastq
    path("*.process.yml"), emit: versions

    script:

    // Base ID name
    sampleName = fastq.getBaseName().replaceAll(~/\.fastq.*$/, '')
    
    // Fastq input can either be a directory or a set of fastq files
    //  Outputs are the same then after allowing a hopefully streamlined pipeline
    if ( fastq.isDirectory() ) {
        """
        artic guppyplex --min-length ${params.min_length} --max-length ${params.max_length} --prefix ${sampleName}_size_selected --directory $fastq

        # Versions #
        cat <<-END_VERSIONS > sizeselect.process.yml
            "${task.process}":
                artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
        END_VERSIONS
        """
    } else {
        """
        mkdir -p input_fastq
        mv $fastq input_fastq/
        artic guppyplex --min-length ${params.min_length} --max-length ${params.max_length} --prefix ${sampleName}_size_selected --directory input_fastq

        # Versions #
        cat <<-END_VERSIONS > sizeselect.process.yml
            "${task.process}":
                artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
        END_VERSIONS
        """
    }
}

process compositeMappingMM2 {
    label 'mediumMem'
    tag { sampleName }

    input:
    tuple val(sampleName), path(singular_fastq), path(composite_ref)

    output:
    tuple val(sampleName), path("${sampleName}.sorted.bam"), emit: comp_bam
    path("*.process.yml"), emit: versions

    script:
    """
    minimap2 -ax map-ont $composite_ref $singular_fastq | samtools view -b | samtools sort -T "temp" -O BAM -o ${sampleName}.sorted.bam

    # Versions #
    cat <<-END_VERSIONS > composite.process.yml
        "${task.process}":
            minimap2: \$(echo \$(minimap2 --version))
            samtools: \$(echo \$(samtools --version | head -n 1 | grep samtools | sed 's/samtools //'))
    END_VERSIONS
    """
}

process removeHumanReads {
    label 'smallCPU'
    tag { sampleName }

    input:
    tuple val(sampleName), path(sorted_bam)

    output:
    tuple val(sampleName), path("${sampleName}.host_removed.sorted.bam"), optional: true, emit: bam
    path "${sampleName}*.csv", emit: csv
    path("*.process.yml"), emit: versions

    script:

    def rev = workflow.commitId ?: workflow.revision ?: workflow.scriptId

    """
    samtools index $sorted_bam
    dehost_nanopore.py --file $sorted_bam --min_reads ${params.min_read_count} --keep_id ${params.keep_ref_id} --output ${sampleName}.host_removed.sorted.bam --revision ${rev}

    # Versions #
    cat <<-END_VERSIONS > dehostedbam.process.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version | head -n 1 | grep samtools | sed 's/samtools //'))
            dehost_nanopore.py: 0.1.0
    END_VERSIONS
    """
}

process regenerateFastqFiles {
    publishDir "${params.outdir}/${params.run_name}/run/fastq_pass/${sampleName}", pattern: "*.host_removed.fastq", mode: "copy"

    label 'smallCPU'
    tag { sampleName }

    input:
    tuple val(sampleName), path(dehosted_bam)

    output:
    tuple val(sampleName), file("${sampleName}.host_removed.fastq"), emit: dehosted_fastq
    path("*.process.yml"), emit: versions

    script:
    """
    samtools fastq $dehosted_bam > ${sampleName}.host_removed.fastq
    # Versions #
    cat <<-END_VERSIONS > dehostedfastq.process.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version | head -n 1 | grep samtools | sed 's/samtools //'))
    END_VERSIONS
    """
}

process regenerateFastqFilesFlat {
    publishDir "${params.outdir}/${params.run_name}/run/fastq_pass/", pattern: "*.host_removed.fastq", mode: "copy"
    label 'smallCPU'
    tag { sampleName }

    input:
    tuple val(sampleName), path(dehosted_bam)

    output:
    tuple val(sampleName), file("${sampleName}.host_removed.fastq"), emit: dehosted_fastq
    path("*.process.yml"), emit: versions

    script:
    """
    samtools fastq $dehosted_bam > ${sampleName}.host_removed.fastq

    # Versions #
    cat <<-END_VERSIONS > dehostedfastq.process.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version | head -n 1 | grep samtools | sed 's/samtools //'))
    END_VERSIONS
    """
}

process regenerateFast5s_MM2 {
    publishDir "${params.outdir}/${params.run_name}/run", pattern: "fast5_pass/${sampleName}", mode: "copy"

    label 'regenerateFast5s'
    tag { sampleName }

    input:
    tuple val(sampleName), path(dehosted_fastq_file), path(fast5_in)

    output:
    path "fast5_pass/${sampleName}", emit: dehosted_fast5
    path("*.process.yml"), emit: versions

    script:
    """
    mkdir -p $sampleName
    mv $dehosted_fastq_file $sampleName/
    bash fast5-dehost-regenerate.sh $sampleName/ $sampleName $fast5_in ${task.cpus}

    # Versions #
    cat <<-END_VERSIONS > regeneratefast5.process.yml
        "${task.process}":
            awk: \$(echo \$(awk --version | head -n 1 | sed 's/GNU Awk //; s/,.*\$//'))
            fast5_subset: NA
    END_VERSIONS
    """
}
