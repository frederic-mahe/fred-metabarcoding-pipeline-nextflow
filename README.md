# fred-metabarcoding-pipeline-nextflow

Fred's metabarcoding pipeline with Nextflow

This is a first attempt at converting to
[Nextflow](https://www.nextflow.io/) the pipeline I use for my own
metabarcoding projects.

Current status:
- [x] draft of the first section of the pipeline (process individual samples),
- [ ] collect info from log files,
- [ ] add checks:
  - [ ] user-set parameters (empty strings, unrealistic values?),
  - [ ] dependencies and minimal versions (vsearch, cutadapt, swarm, bash >= 4),
- [ ] adapt to work with slurm,
- [ ] automatically deduce the fastq file naming pattern,
- [ ] automatically deduce compression (gz, bz2) or the lack-of,
- [ ] extend to multiplexed datasets,
- [ ] draft of the second section of the pipeline (work at the study scale)

## Run workflow with test dataset

```bash
nextflow run -r main frederic-mahe/fred-metabarcoding-pipeline-nextflow -profile test
```