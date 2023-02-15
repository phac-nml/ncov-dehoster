## Version 0.3.0
----------------
General Developer Changes:
- Changed how conda is implemented to match nextflow version updates
- Fixed tests

## Version 0.2.0
----------------
Argument Changes:
- Added `--keep_ref_id <ID>` as an argument for illumina and nanopore fastq data
    - Allows for a different reference ID to be kept if the reference sequence is changed
- Removed `--covid_ref_id <ID>` as an illumina specific argument
    - Replaced by above
- Removed the `illumina_threads` and `nanopore_threads` resource arguments that were not super useful
    - Changed to using `task.cpus` to keep everything in line
    - Resources should be changed using the process name or the process tag now

Parameter Changes:
- Default minimum read length for nanopore data (`--min_length`) set to 350 from 400 to better accommodate shorter amplicon schemes
- Help command updated

Output Changes:
- New output file: `process_versions.yml`
    - File contains process name along with tools and their versions used in it

General Developer Changes:
- Changed code order for nanopore workflows
    - All minimap2 code should be above nanostripper in the files now
- General reformatting of code to try to make everything more consistent/cohesive
- Moved the minimap2 workflow above the nanostripper workflow in all processes

## Version 0.1.0
----------------
Argument Changes:
- Nanopore Nanostripper Pipeline:
    - Accessed with `--nanopore` and `--nanostripper` flags
    - Fast5 input changed from `--directory` to `--fast5_directory`

Added Functionalities:
- Added nanopore fastq host removal pipeline focused on minimap2
    - dehost_nanopore.py to accompany this
- Added help statement module that is printed when `--help` is passed 
- Added pysam and pip to nanopore environment

General Developer Changes:
- DSL2 enabled (no longer on DSL2 preview)
    - Small changes to some modules to address this but no functionality changes
    - Nextflow version >= 21.04.0 needed
- Adjusted labels for CPU and MEM usage
    - mediumMem replaces mediumCPU
- Adjusted config file structure
- Adjusted how modules are imported
- Adjusted some nanopore module names to not overlap
    - nanoporeDehosting  -> nanoporeNanostripperDehosting
    - fastqSizeSelection -> fastqSizeSelection_NS 
    - regenerateFast5s   -> regenerateFast5s_NS


## Version 0.0.1
----------------

- Initial Release that works
