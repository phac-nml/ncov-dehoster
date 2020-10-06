# ncov-dehoster

## About:
Nextflow pipeline that removes human reads from SARS-CoV-2 Illumina and Nanopore (not ready) sequencing data.

## Run:
Run the Illumina pipeline with the simple command:

```
nextflow run phac-nml/ncov-dehoster -profile conda --directory <path/to/paired_reads/dir> --illumina
```

If you want to modify the executor (say use slurm for example) or the resource allocations generate and
add a custom profile that you can specify with `-profile`

## Notes:

- If you have the dependencies installed already you can skip the `-profile conda` step

- This is an initial version and there are potentially still bugs and changes to be made
