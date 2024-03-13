process seqtkRandomSubsample {
    // Output logic for illumina and nanopore flat/normal structures
    publishDir = [
        path: { (params.illumina ? "${params.outdir}/dehosted_paired_fastqs" : params.flat ? "${params.outdir}/${params.run_name}/run/fastq_pass/" : "${params.outdir}/${params.run_name}/run/fastq_pass/${sampleName}") },
        pattern: "*.fastq*",
        mode: "copy"
    ]
    label 'process_low'
    tag { sampleName }

    input:
    tuple val(sampleName), path(reads)
    val(sample_size)

    output:
    tuple val(sampleName), path("*.fastq*"), emit: reads
    path "versions.yml", emit: versions

    script:
    if ( params.illumina ) {
        MID_EXT = "_downsampled"
        GZIP_ARG = "| gzip --no-name"
    } else {
        MID_EXT = ".downsampled"
        GZIP_ARG = ""
    }
    """
    # For loop to capture both nanopore and illumina reads
    for f in $reads;
    do
        FINAL_EXT=\$(echo \$f | sed 's/${sampleName}//')
        seqtk \\
            sample \\
            -s ${params.downsample_seed} \\
            \$f \\
            $sample_size \\
        $GZIP_ARG \\
        > ${sampleName}${MID_EXT}\$FINAL_EXT
    done

    # Versions #
    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
process samtoolsAmpliconDownsample {
    publishDir "${params.outdir}/downsample_stats", pattern: "${sampleName}_region_counts.csv", mode: "copy"
    label 'process_low'
    tag { sampleName }

    input:
    tuple val(sampleName), path(bam)
    path amplicons
    val platform
    val sample_size
    val seed

    output:
    tuple val(sampleName), path("${sampleName}_downsampled.bam"), emit: bam
    tuple val(sampleName), path("${sampleName}_region_counts.csv"), emit: csv
    path "versions.yml", emit: versions

    script:
    """
    # Get approximate per amplicon read count
    perAmpliconReadCount=\$((${sample_size}/\$(wc -l ${amplicons} | cut -f1 -d' ')))

    # Files must be sorted to work with downsampling
    samtools sort ${bam} > ${sampleName}.sorted.bam
    samtools index ${sampleName}.sorted.bam

    # Downsample with samtools view and awk
    bash downsample_bam.sh \\
        --bam ${sampleName}.sorted.bam \\
        -a ${amplicons} \\
        -p ${platform} \\
        --read-count \$perAmpliconReadCount \\
        --seed ${seed} \\
        --sample-name ${sampleName}

    # Versions #
    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version | head -n 1 | grep samtools | sed 's/samtools //'))
    END_VERSIONS
    """
}
