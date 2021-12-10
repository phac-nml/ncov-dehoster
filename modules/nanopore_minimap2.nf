process generateMinimap2Index {
    publishDir "${params.outdir}/compositeMinimapIndex", pattern: "composite_ref.mmi", mode: "symlink"

    label 'bwa_mem'

    input:
    path(human_ref)
    path(viral_ref)

    output:
    file "composite_ref.mmi"

    script:
    """
    cat $human_ref $viral_ref > composite_reference.fa
    minimap2 -d composite_ref.mmi composite_reference.fa
    """
}

process guppyplexSizeSelection {
    publishDir "${params.outdir}/fastqSizeSelect", pattern: "*.fastq", mode: "copy"

    label 'smallcpu'

    input:
    path(fastq)

    output:
    tuple val(sampleName), path("${sampleName}_size_selected*.fastq")

    script:

    // Base ID name
    sampleName = fastq.getBaseName()
    
    // Fastq input can either be a directory or a set of fastq files
    //  Outputs are the same then after
    if ( fastq.isDirectory() ) {
        """
        artic guppyplex --min-length ${params.min_length} --max-length ${params.max_length} --prefix ${sampleName}_size_selected --directory $fastq
        """
    } else {
        """
        mkdir -p input_fastq
        mv $fastq input_fastq/
        artic guppyplex --min-length ${params.min_length} --max-length ${params.max_length} --prefix ${sampleName}_size_selected --directory input_fastq
        """
    }
}

process compositeMapping {
    publishDir "${params.outdir}/compositeMapping", pattern: "*.sorted.bam", mode: "symlink"

    label 'smallMem'

    input:
    tuple val(sampleName), path(singular_fastq), path(composite_ref)

    output:
    tuple val(sampleName), path("${sampleName}.sorted.bam")

    script:
    """
    minimap2 -ax map-ont $composite_ref $singular_fastq | samtools view -b | samtools sort -T "temp" -O BAM -o ${sampleName}.sorted.bam
    """
}

process removeHumanReads {
    publishDir "${params.outdir}/host_removed", pattern: "*.sorted.bam", mode: "symlink"

    label 'smallMem'

    input:
    tuple val(sampleName), path(sorted_bam)

    output:
    tuple val(sampleName), path("${sampleName}.host_removed.sorted.bam")

    script:
    """
    samtools index $sorted_bam
    dehost_nanopore.py --file $sorted_bam --output ${sampleName}.host_removed.sorted.bam
    """
}

process generateFastqFiles {
    publishDir "${params.outdir}/final_fastq", pattern: "*.host_removed.fastq", mode: "copy"

    label 'smallMem'

    input:
    tuple val(sampleName), path(dehosted_bam)

    output:
    tuple val(sampleName), file("${sampleName}.host_removed.fastq")

    script:
    """
    samtools fastq $dehosted_bam > ${sampleName}.host_removed.fastq
    """
}

process regenerateDehostedFast5s {
    publishDir "${params.outdir}/", pattern: "fast5_pass/${sampleName}", mode: "copy"

    input:
    tuple val(sampleName), path(dehosted_fastq_file), path(fast5_in)

    output:
    path "fast5_pass/${sampleName}"

    script:
    """
    mkdir -p $sampleName
    mv $dehosted_fastq_file $sampleName/
    bash fast5-dehost-regenerate.sh $sampleName/ $sampleName $fast5_in
    """
}

// process regenerateDehostedSequencingSummary {

// }
