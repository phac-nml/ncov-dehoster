process nanostripper {

    publishDir "${params.outdir}/${params.run_name}/run", pattern: "fast5_dehosted/${barcodeName}", mode: "copy"

    label 'nanostripper'

    tag { barcodeName }

    input:
    tuple path(folder), file(sars_reference), file(human_reference)

    output:
    path "fast5_dehosted/${barcodeName}", emit: dehostedFast5
    val "${barcodeName}", emit: barcode

    script:

    barcodeName = folder.getBaseName().replaceAll(~/\.*$/, '')

    // Temporary path to nanostripper used while conda env in progress
    // Temporary path to nanostripper tool while I look for a better solution
    """
    ${params.nanostripper_tool_path}nanostripper -out ./fast5_dehosted -t 10 ${sars_reference} ${human_reference} ${folder}
    mv ./fast5_dehosted/nanostripper_summary.txt ./fast5_dehosted/${barcodeName}
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

    publishDir "${params.outdir}/${params.run_name}/run", pattern: "fastq_pass", mode: "copy"

    input:
    file(dehosted_combined_fastq)

    output:
    path "fastq_pass"
    path "fastq_pass/*", type: 'dir', emit: barcodes

    script:
    """
    guppy_barcoder --require_barcodes_both_ends -i ./ -s fastq_pass --arrangements_files 'barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg barcode_arrs_nb96.cfg'
    """
}

process combineFast5Barcodes {

    // Done to make a single directory of the fast5 dehosted only files for regenerateFast5 process

    label 'smallmem'

    input:
    path(dehosted_fast5s)

    output:
    path "fast5_pass_dehosted_only"

    script:

    """
    mkdir -p fast5_pass_dehosted_only
    mv $dehosted_fast5s fast5_pass_dehosted_only/
    """
}

process regenerateFast5s {

    publishDir "${params.outdir}/${params.run_name}/run", pattern: "fast5_pass/${barcodeName}", mode: "copy"

    input:
    path(dehosted_fastq_barcode)
    path(fast5_dehosted)

    output:
    path "fast5_pass/${barcodeName}" , emit: fast5_pass
    path "${barcodeName}.csv" , emit: csv

    script:

    barcodeName = dehosted_fastq_barcode.getBaseName().replaceAll(~/\.*$/, '')

    def rev = workflow.commitId ?: workflow.revision ?: workflow.scriptId

    """
    bash fast5-dehost-regenerate.sh $dehosted_fastq_barcode $barcodeName $fast5_dehosted

    bash generate-csv.sh ${barcodeName} ${fast5_dehosted}/${barcodeName}/nanostripper_summary.txt fast5_pass/${barcodeName}/filename_mapping.txt ${rev}
    """
}

process generateSimpleSequencingSummary {

    publishDir "${params.outdir}/${params.run_name}/run", pattern: "sequencing_summary.txt", mode: "copy"

    input:
    path(fast5_barcode)

    output:
    path "sequencing_summary.txt"

    script:
    """
    cat <(echo -e "read_id\tfilename") <(cat ./*/filename_mapping.txt) > sequencing_summary.txt
    """
}

process combineCSVs {

    publishDir "${params.outdir}/${params.run_name}", pattern: "removal_summary.csv", mode: "copy"

    label 'smallcpu'

    input:
    path(csvs)

    output:
    path("removal_summary.csv")

    script:
    """
    csvtk concat *.csv > summary.csv
    csvtk sort -k1 summary.csv > removal_summary.csv
    """
}
