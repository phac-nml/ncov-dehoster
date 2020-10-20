process nanostripper {

    publishDir "${params.outdir}/${params.run_name}", pattern: "fast5_dehosted/${barcodeName}", mode: "copy"

    label 'nanostripper'

    tag { barcodeName }

    input:
    tuple path(barcode), file(sars_reference), file(human_reference)

    output:
    tuple barcodeName, path("fast5_dehosted/${barcodeName}")

    script:

    barcodeName = barcode.getBaseName().replaceAll(~/\.*$/, '')

    """
    nanostripper -out ./fast5_dehosted -t 10 ${sars_reference} ${human_reference} ${barcode} 
    """
}

process guppyBasecaller {

    publishDir "${params.outdir}/${params.run_name}/test", pattern: "*.txt", mode: "copy"

    label 'guppy'

    // input:
    // tuple path(barcode), file(sars_reference), file(human_reference)

    // output:
    // tuple sampleName, file("${sampleName}.sorted.bam")

    script:

    """
    echo hi >> f.txt
    """
}
