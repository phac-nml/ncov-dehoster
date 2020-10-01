process zeroReads {

    label 'smallcpu'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "out.txt", mode: "copy"

    input:
    file(fastq)

    output:
    file('out.txt')
    

    script:
    """
    echo ${fastq} >> out.txt
    """
}