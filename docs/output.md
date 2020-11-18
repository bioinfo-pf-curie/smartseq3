# Outputs

<!-- TODO update with the output of your pipeline -->

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes the data using the steps presented in the main README file.  
Briefly, its goal is to process <!-- TODO --> data for any protocol, with or without control samples, and with or without spike-ins.

The directories listed below will be created in the output directory after the pipeline has finished. 

## Sequencing quality

### FastQC
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

**Output directory: `fastqc`**

* `sample_fastqc.html`
  * FastQC report, containing quality metrics for your untrimmed raw fastq files.
* `zips/sample_fastqc.zip`
  * zip file containing the FastQC report, tab-delimited data file and plot images.

## Read mapping

### Alignment

Different tools can be used for read alignment (`STAR`, `BWA-mem`, `Bowtie2`). The mapping statistics (`Total Reads`, `Aligned Reads`, `Unique Reads`, `Multiple Reads`) are also presented in the main summary table.

> **NB:** by default, one alignment is randomly reported in case of multiple mapping sites. If necessary, these reads can be filtered using the `--mapq` option. In addition, in case of paired-end sequencing reads, singleton are discarded from the analysis.

**Output directory: `mapping`**

* `sample.bam`
  * Aligned reads. Available only if (`--saveAlignedIntermediates`) is used.
* `sample_sorted.bam`
  * Aligned reads sorted by chromosome position. Available only if (`--saveAlignedIntermediates`) is used.
* `sample_sorted.bam.bai`
  * Index of aligned and sorted reads. Available only if (`--saveAlignedIntermediates`) is used.

The mapping statistics are presented in the MultiQC report as follows.  
In general, we expect more than 80% of aligned reads. Samples with less than 50% of mapped reads should be further investigated, and checked for adapter content, contamination, etc.

![MultiQC - Bowtie2 stats plot](images/bowtie2.png)

### Duplicates

[Picard MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html) is used to mark and remove the duplicates. 
The results are presented in the `General Metrics` table. Duplicate reads are **removed** by default from the aligned reads to mitigate for fragments in the library that may have been sequenced more than once due to PCR biases. There is an option to keep duplicate reads with the `--keepDups` parameter but it is generally recommended to remove them to avoid the wrong interpretation of the results.	

**Output directory: `mapping`**

* `sample_filtered.bam`
  * Aligned reads after filtering (`--mapq`, `--keepDups`).
* `sample_filtered.bam.bai`
  * Index of aligned reads after filtering.

![MultiQC - Picard MarkDup stats plot](images/picard_deduplication.png)

## Quality controls

From the filtered and aligned read files, the pipeline runs several quality control steps presented below.

### Sequencing complexity

The [Preseq](http://smithlabresearch.org/software/preseq/) package aims at predicting and estimating the complexity of a genomic sequencing library, equivalent to predicting and estimating the number of redundant reads from a given sequencing depth and how many will be expected from additional sequencing using an initial sequencing experiment. The estimates can then be used to examine the utility of further sequencing, optimize the sequencing depth, or to screen multiple libraries to avoid low complexity samples. The dashed line shows a perfectly complex library where total reads = unique reads. Note that these are predictive numbers only, not absolute. The MultiQC plot can sometimes give extreme sequencing depth on the X axis - click and drag from the left side of the plot to zoom in on more realistic numbers.

**Output directory: `preseq`**

* `sample_ccurve.txt`
  * Preseq expected future yield file.

![MultiQC - Preseq library complexity plot](images/preseq_plot.png)

## MultiQC
[MultiQC](http://multiqc.info) is a visualisation tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available within the report data directory.

The pipeline has special steps which allow the software versions used to be reported in the MultiQC output for future traceability.

**Output directory: `multiqc`**

* `multiqc_report.html`
  * MultiQC report - a standalone HTML file that can be viewed in your web browser.
* `multiqc_data/`
  * Directory containing parsed statistics from the different tools used in the pipeline.

For more information about how to use MultiQC reports, see http://multiqc.info.
