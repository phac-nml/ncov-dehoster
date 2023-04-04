process seqtkSubsample {
    publishDir = [
        path: { (params.illumina ? "${params.outdir}/dehosted_paired_fastqs" : params.flat ? "${params.outdir}/${params.run_name}/run/fastq_pass/" : "${params.outdir}/${params.run_name}/run/fastq_pass/${sampleName}") },
        pattern: "*.fastq*",
        mode: "copy"
    ]

    tag { sampleName }
    label 'smallCPU'

    conda "bioconda::seqtk=1.3"

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
