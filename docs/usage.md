# Usage

## Table of contents

* [Introduction](#general-nextflow-info)
* [Running the pipeline](#running-the-pipeline)
* [Main arguments](#main-arguments)
    * [`--reads`](#--reads)
    * [`--samplePlan`](#--samplePlan)
* [Reference genomes](#reference-genomes)
    * [`--genome`](#--genome)
    * [`--genomeAnnotationPath`](#--genomeAnnotationPath)
* [Nextflow profiles](#nextflow-profiles)
* [Job resources](#job-resources)
* [Other command line parameters](#other-command-line-parameters)
    * [`--skip*`](#--skip*)
    * [`--metadata`](#--metadta)
    * [`--outDir`](#--outdir)
    * [`--email`](#--email)
    * [`-name`](#-name-single-dash)
    * [`-resume`](#-resume-single-dash)
    * [`-c`](#-c-single-dash)
    * [`--maxMemory`](#--maxMemory)
    * [`--maxTime`](#--maxTime)
    * [`--maxCpus`](#--maxCpus)
    * [`--multiqcConfig`](#--multiqcConfig)

## General Nextflow info
Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline
The typical command for running the pipeline is as follows:
```bash
nextflow run main.nf --reads '*_R{1,2}.fastq.gz' -profile 'singularity'
```

This will launch the pipeline with the `singularity` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

You can change the output director using the `--outDir/-w` options.

## Main arguments

### `--reads`
Use this to specify the location of your input FastQ files. For example:

```bash
--reads 'path/to/data/sample_*_{1,2}.fastq'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs.

If left unspecified, a default pattern is used: `data/*{1,2}.fastq.gz`


### `--samplePlan`
Use this to specify a sample plan file instead of a regular expression to find fastq files. For example :

```bash
--samplePlan 'path/to/data/sample_plan.csv
```

The sample plan is a csv file with the following information :

Sample ID | Sample Name | Path to R1 fastq file | Path to R2 fastq file

## Reference Genomes

All information about genomes and annotation are available in `conf/genomes.config`.

### `-genome`

There are different species supported in the genomes references file. To run the pipeline, you must specify which to use with the `--genome` flag.

You can find the keys to specify the genomes in the [genomes config file](../conf/genomes.config). Common genomes that are supported are:

* Human
  * `--genome hg38`
  * `--genome hg19`
* Mouse
  * `--genome mm10`
  * `--genome mm9`

> There are numerous others - check the config file for more.

Note that you can use the same configuration setup to save sets of reference files for your own use, even if they are not part of the genomes resource. 
The syntax for this reference configuration is as follows:

```nextflow
params {
  genomes {
    'hg19' {
      starIndex = '<path to the STAR index files>'
      gtf       = '<path to a gene annotation file in gtf format>'
	    bed12     = '<path to annotation file in b bed12 format>'
    }
    // Any number of additional genomes, key is used with --genome
  }
}
```

Note that these paths can be updated on command line using the following parameters:
- `--starIndex` - Path to STAR index
- `--gtf` - Path to GTF file
- `--bed12` - Path to bed12 file

### `--genomeAnnotationPath`

The `--genomeAnnotationPath` define where all annotations are stored. This path can be defined on the command line or set up in the different configuration file during the pipeline installation.
See `conf/installation.md` for details.

## Nextflow profiles

Different Nextflow profiles can be used. See [Profiles](profiles.md) for details.

## Job resources

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time (see the [`conf/process.conf`](../conf/process.config) file). 
For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If 
it still fails after three times then the pipeline is stopped.

## Other command line parameters

### `--skip*`

The pipeline is made with a few *skip* options that allow to skip optional steps in the workflow.
The following options can be used:
- `--skipMultiqc` - Skip MultiQC
- - `--skipSoftVersion` - Skip showing software version in MultiQC
- `--skipGeneCov` - Skip gene body coverage

### `--metadata`

Specify a two-columns (tab-delimited) metadata file to diplay in the final Multiqc report.

### `--outdir`

The output directory where the results will be saved.

### `--email`

Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to speicfy this on the command line for every run.

### `-name`

Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`

Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

### `--maxMemory`

Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--maxMemory '8.GB'`

### `--maxTime`

Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--maxTime '2.h'`

### `--maxCpus`

Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--maxCpus 1`

### `--multiqcConfig`

Specify a path to a custom MultiQC configuration file.

## Job resources

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time (see the `conf/base.conf` file). 
For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.
