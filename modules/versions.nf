process outputVersionsIllumina {
    publishDir "${params.outdir}", pattern: "process_versions.yml", mode: "copy"

    label 'smallCPU'

    input:
    path versions

    output:
    path "process_versions.yml"

    script:
    def rev = workflow.commitId ?: workflow.revision ?: workflow.scriptId
    """
    echo '"ncov-dehoster-${rev}":' > process_versions.yml
    cat *.process.yml >> process_versions.yml
    """
}
process outputVersionsNanopore {
    publishDir "${params.outdir}/${params.run_name}", pattern: "process_versions.yml", mode: "copy"

    label 'smallCPU'

    input:
    path versions

    output:
    path "process_versions.yml"

    script:
    def rev = workflow.commitId ?: workflow.revision ?: workflow.scriptId
    """
    echo '"ncov-dehoster-${rev}":' > process_versions.yml
    cat *.process.yml >> process_versions.yml
    """
}
