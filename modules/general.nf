process seqtkSubsample {
    if ( params.illumina ) {
        publishDir "${params.outdir}/dehosted_paired_fastqs", pattern: "*.fastq.gz", mode: "copy"
    } else if ( params.flat ) {
        publishDir "${params.outdir}/${params.run_name}/run/fastq_pass/", pattern: "*.fastq.gz", mode: "copy"
    } else {
        publishDir "${params.outdir}/${params.run_name}/run/fastq_pass/${sampleName}", pattern: "*.fastq.gz", mode: "copy"
    }

    tag { sampleName }
    label 'smallCPU'

    conda "bioconda::seqtk=1.3"

    input:
    tuple val(sampleName), path(reads)
    val(sample_size)

    output:
    tuple val(sampleName), path("*.fastq.gz"), emit: reads
    path "versions.yml", emit: versions

    script:
    if ( params.illumina ) {
        MID_EXT = "_downsampled"
    } else {
        MID_EXT = ".downsampled"
    }
    """
    for f in $reads;
    do
        FINAL_EXT=\$(echo \$f | sed 's/${sampleName}//')
        seqtk \
            sample \
            -s ${params.downsample_seed} \
            \$f \
            $sample_size \
        | gzip --no-name \
        > ${sampleName}${MID_EXT}\$FINAL_EXT
    done

    # Versions #
    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
