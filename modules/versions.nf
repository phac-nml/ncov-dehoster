process outputVersions {
    publishDir "${params.outdir}", pattern: "process_versions.yml", mode: "copy"
    label 'process_low'

    input:
    path versions

    output:
    path "process_versions.yml"

    script:
    def rev = workflow.commitId ?: workflow.revision ?: workflow.scriptId
    """
    echo '"ncov-dehoster-${rev}":' > process_versions.yml
    cat $versions >> process_versions.yml
    """
}
