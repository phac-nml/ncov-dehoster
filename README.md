# ncov-dehoster

## About:
Nextflow pipeline that removes human reads from SARS-CoV-2 Illumina and Nanopore (not ready) sequencing data.

Illumina - Takes a competitive mapping approach with bwa mem to remove human reads from the input fastq files

## Run:
Run the Illumina pipeline with the simple command:

```
nextflow run phac-nml/ncov-dehoster -profile conda --directory <path/to/paired_reads/dir> --illumina --human_ref <path/to/reference>
```

If you want to modify the executor (say use slurm for example) or the resource allocations generate and
add a custom profile that you can specify with `-profile`

## Notes:

- If you have the dependencies installed already you can skip the `-profile conda` step

- You will need a human reference genome and to specify it with `--human_ref <path/to/reference>`. The nml path is in by default (sorry)

- Slowest part of the pipeline is indexing the composite human/Sars CoV2 reference with bwa index, if you have a bwa index of this
pass the folder containing the indexes to the script with the parameter `--composite_bwa_index <path/to/bwa_indexes/>`

- This is an initial version and there are potentially still bugs and changes to be made

## Upcoming (Hopefully) Additions and Changes:

- Nanopore removal and re-basecalling of fast5 raw data

- More info in `removal_summary.csv` output
