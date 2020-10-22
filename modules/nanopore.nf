process nanostripper {

    publishDir "${params.outdir}/${params.run_name}", pattern: "fast5_dehosted/${barcodeName}", mode: "copy"

    label 'nanostripper'

    tag { barcodeName }

    input:
    tuple path(folder), file(sars_reference), file(human_reference)

    output:
    path "fast5_dehosted/${barcodeName}", emit: dehostedFast5
    val "${barcodeName}", emit: barcode

    script:

    barcodeName = folder.getBaseName().replaceAll(~/\.*$/, '')

    // Temporary path to nanostripper while conda env in progress
    """
    /Drives/W/Projects/covid-19/nml_reads/vdd/dehosting_runs/nanostripper/nanostripper -out ./fast5_dehosted -t 10 ${sars_reference} ${human_reference} ${folder} 
    """
}

process guppyBasecallerGPU {

    label 'guppyGPU'

    input:
    path(dehosted_fast5s)

    output:
    path "fastq_pass_dehosted_only"

    script:

    """
    bash guppy-gpu.sh ${params.max_parallel_basecalling} ${params.gpu_per_node}
    """
}

process guppyBasecallerCPU {

    label 'guppyCPU'

    input:
    path(dehosted_fast5_barcode)

    output:
    path "fastq_pass_dehosted_only/$dehosted_fast5_barcode"

    script:

    """
    guppy_basecaller -c dna_r9.4.1_450bps_hac.cfg -r -i $dehosted_fast5_barcode -s fastq_pass_dehosted_only/$dehosted_fast5_barcode
    """
}

process combineFastq {

    label 'largeMem'

    input:
    path(fastq_pass_dehosted_only)

    output:
    file "combined.fastq"

    script:
    """
    cat ./${fastq_pass_dehosted_only}/*/*.fastq > combined.fastq
    """
}

process fastqSizeSelection {

    label 'largeMem'

    input:
    file(combined_fastq)

    output:
    file "${params.run_name}*.fastq"

    script:
    """
    mkdir -p fastq_pass_combined_dehosted
    mv $combined_fastq fastq_pass_combined_dehosted
    artic guppyplex --min-length ${params.min_length} --max-length ${params.max_length} --directory ./fastq_pass_combined_dehosted --prefix ${params.run_name} 
    """
}

process fastqDemultiplex {

    label 'largeMem'

    publishDir "${params.outdir}/${params.run_name}", pattern: "fastq_pass", mode: "copy"

    input:
    file(dehosted_combined_fastq)

    output:
    path "fastq_pass"

    script:
    """
    guppy_barcoder --require_barcodes_both_ends -i ./ -s fastq_pass --arrangements_files 'barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg barcode_arrs_nb96.cfg'
    """
}
