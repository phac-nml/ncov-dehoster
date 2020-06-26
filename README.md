# ncov-dehumanizer

## About:
Nextflow and slurm pipeline that removes human reads from SARS-CoV-2 nanopore sequencing data

## Run:
Run the pipeline with the simple command:

```
nextflow run phac-nml/ncov-dehumanizer -profile conda --directory <path/to/reads/dir>
```

## Notes:

- If you have the dependencies installed already you can skip the `-profile conda` step

- Set up to only really work on the NML cluster for the moment

- This is an initial version and there are potentially still bugs and changes to be made
