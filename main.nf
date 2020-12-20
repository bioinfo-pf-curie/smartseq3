#!/usr/bin/env nextflow

/*
Copyright Institut Curie 2020
This software is a computer program whose purpose is to analyze high-throughput sequencing data.
You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND.
Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data.
The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.

This script is based on the nf-core guidelines. See https://nf-co.re/ for more information
*/

/*
========================================================================================
SmartSeq3
========================================================================================
 #### Homepage / Documentation
https://gitlab.curie.fr/sc-platform/smartseq3
----------------------------------------------------------------------------------------
*/

// File with text to display when a developement version is used
devMessageFile = file("$baseDir/assets/devMessage.txt")

def helpMessage() {
  if ("${workflow.manifest.version}" =~ /dev/ ){
    devMess = file("$baseDir/assets/devMessage.txt")
    log.info devMessageFile.text
  }

  log.info """
  SmartSeq3 v${workflow.manifest.version}
  ======================================================================

  Usage:
  nextflow run main.nf --reads '*_R{1,2}.fastq.gz' --genome 'hg19' -profile conda
  nextflow run main.nf --samplePlan samplePlan --genome 'hg19' -profile conda

  Mandatory arguments:
    --reads [file]                Path to input data (must be surrounded with quotes)
    --samplePlan [file]           Path to sample plan input file (cannot be used with --reads)
    --genome [str]                Name of genome reference
    -profile [str]                Configuration profile to use. test / conda / multiconda / path / multipath / singularity / docker / cluster (see below)
  
  Inputs:
    --design [file]               Path to design file for extended analysis  
    --singleEnd [bool]            Specifies that the input is single-end reads

  Options:
      --minCountPerCell1             First minimum umi counts per cell. Default: 500
      --minCountPerCell2             Second minimum umi counts per cell. Default: 1000
      --minCountPerCellGene1         First minimum gene counts per cell. Default: 100
      --minCountPerCellGene2         Second minimum gene counts per cell. Default: 200

  Skip options: All are false by default
    --skipSoftVersion [bool]      Do not report software version
    --skipMultiQC [bool]          Skips MultiQC
  
  Genomes: If not specified in the configuration file or if you wish to overwrite any of the references given by the --genome field
  --genomeAnnotationPath [file]      Path  to genome annotation folder
  --starIndex [dir]                  Index for STAR aligner

  Other options:
    --outdir [file]               The output directory where the results will be saved
    -name [str]                   Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
 
  =======================================================
  Available Profiles

    -profile test                Set up the test dataset
    -profile conda               Build a single conda for with all tools used by the different processes before running the pipeline
    -profile multiconda          Build a new conda environment for each tools used by the different processes before running the pipeline
    -profile path                Use the path defined in the configuration for all tools
    -profile multipath           Use the paths defined in the configuration for each tool
    -profile docker              Use the Docker containers for each process
    -profile singularity         Use the singularity images for each process
    -profile cluster             Run the workflow on the cluster, instead of locally

  """.stripIndent()
}

/**********************************
 * SET UP CONFIGURATION VARIABLES *
 **********************************/

// Show help message
if (params.help){
  helpMessage()
  exit 0
}

// Configurable reference genomes

// TODO - ADD HERE ANY ANNOTATION

if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
  exit 1, "The provided genome '${params.genome}' is not available in the genomes.config file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

//Stage config files
Channel
  .fromPath(params.multiqcConfig, checkIfExists: true)
  .set{chMultiqcConfig}
chOutputDocs = file("$baseDir/docs/output.md", checkIfExists: true)
chOutputDocsImages = file("$baseDir/docs/images/", checkIfExists: true)

//Has the run name been specified by the user?
//This has the bonus effect of catching both -name and --name
customRunName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  customRunName = workflow.runName
}

/************
 * CHANNELS *
 ************/

// Validate inputs
if ((params.reads && params.samplePlan) || (params.readPaths && params.samplePlan)){
  exit 1, "Input reads must be defined using either '--reads' or '--samplePlan' parameters. Please choose one way."
}

if ( params.metadata ){
  Channel
    .fromPath( params.metadata )
    .ifEmpty { exit 1, "Metadata file not found: ${params.metadata}" }
    .set { chMetadata }
}else{
  chMetadata=Channel.empty()
}                 


/* Genomes -------------*/

// Configurable reference genomes
genomeRef = params.genome

params.starIndex = genomeRef ? params.genomes[ genomeRef ].starIndex ?: false : false
if (params.starIndex){
  Channel
    .fromPath(params.starIndex, checkIfExists: true)
    .ifEmpty {exit 1, "STAR index file not found: ${params.starIndex}"}
    .set { chStar }
} else {
  exit 1, "STAR index file not found: ${params.starIndex}"
}

params.gtf = genomeRef ? params.genomes[ genomeRef ].gtf ?: false : false
if (params.gtf) {
  Channel
    .fromPath(params.gtf, checkIfExists: true)
    .into { chGtfSTAR; chGtfFC }
}else {
  exit 1, "GTF annotation file not not found: ${params.gtf}"
}

params.bed12 = genomeRef ? params.genomes[ genomeRef ].bed12 ?: false : false
if (params.bed12) {
  Channel  
    .fromPath(params.bed12)
    .ifEmpty { exit 1, "BED12 annotation file not found: ${params.bed12}" }
    .set { chBedGeneCov } 
}else {
  exit 1, "GTF annotation file not not found: ${params.bed12}"
}

/*----------------*/

// Create a channel for input read files
if(params.samplePlan){
  if(params.singleEnd){
    Channel
      .from(file("${params.samplePlan}"))
      .splitCsv(header: false)
      .map{ row -> [ row[0], [file(row[2])]] }
      .into { rawReadsFastqc; chMergeReadsFastq }
  }else{
    Channel
      .from(file("${params.samplePlan}"))
      .splitCsv(header: false)
      .map{ row -> [ row[0], [file(row[2]), file(row[3])]] }
      .into { rawReadsFastqc ; chMergeReadsFastq}
   }
  params.reads=false
}
else if(params.readPaths){
  if(params.singleEnd){
    Channel
      .from(params.readPaths)
      .map { row -> [ row[0], [file(row[1][0])]] }
      .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied." }
      .into { rawReadsFastqc ; chMergeReadsFastq}
  } else {
    Channel
      .from(params.readPaths)
      .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
      .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied." }
      .into { rawReadsFastqc ; chMergeReadsFastq}
  }
} else {
  Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .into { rawReadsFastqc ; chMergeReadsFastq}
}

// Make sample plan if not available
/**********************************/

if (params.samplePlan){
  Channel
    .fromPath(params.samplePlan)
    .into {chSplan; chSplanCheck}
}else if(params.readPaths){
  if (params.singleEnd){
    Channel
      .from(params.readPaths)
      .collectFile() {
        item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + '\n']
       }
      .into{ chSplan; chSplanCheck }
  }else{
     Channel
       .from(params.readPaths)
       .collectFile() {
         item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + ',' + item[1][1] + '\n']
        }
       .into{ chSplan; chSplanCheck }
  }
} else if(params.bamPaths){
  Channel
    .from(params.bamPaths)
    .collectFile() {
      item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + '\n']
     }
    .into{ chSplan; chSplanCheck }
  params.aligner = false
} else {
  if (params.singleEnd){
    Channel
      .fromFilePairs( params.reads, size: 1 )
      .collectFile() {
         item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + '\n']
      }     
      .into { chSplan; chSplanCheck }
  }else{
    Channel
      .fromFilePairs( params.reads, size: 2 )
      .collectFile() {
         item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + ',' + item[1][1] + '\n']
      }     
      .into { chSplan; chSplanCheck }
   }
}


/***************
 * Design file *
 ***************/

// TODO - UPDATE BASED ON YOUR DESIGN

if (params.design){
  Channel
    .fromPath(params.design)
    .ifEmpty { exit 1, "Design file not found: ${params.design}" }
    .into { chDesignCheck; chDesignControl; chDesignMqc }

  chDesignControl
    .splitCsv(header:true)
    .map { row ->
      if(row.CONTROLID==""){row.CONTROLID='NO_INPUT'}
      return [ row.SAMPLEID, row.CONTROLID, row.SAMPLENAME, row.GROUP, row.PEAKTYPE ]
     }
    .set { chDesignControl }

  // Create special channel to deal with no input cases
  Channel
    .from( ["NO_INPUT", ["NO_FILE","NO_FILE"]] )
    .toList()
    .set{ chNoInput }
}else{
  chDesignControl = Channel.empty()
  chDesignCheck = Channel.empty()
  chDesignMqc = Channel.empty()
}


/*******************
 * Header log info *
 *******************/

if ("${workflow.manifest.version}" =~ /dev/ ){
   log.info devMessageFile.text
}

log.info """=======================================================

smartSeq3 v${workflow.manifest.version}
======================================================="""
def summary = [:]

/* ADDED -------------*/
summary['Pipeline Name']  = 'SmartSeq3'
summary['Pipeline Version'] = workflow.manifest.version
//summary['Run Name']     = custom_runName ?: workflow.runName
summary['Command Line'] = workflow.commandLine
if (params.samplePlan) {
   summary['SamplePlan']   = params.samplePlan
}else{
   summary['Reads']        = params.reads
}
summary['Genome']       = params.genome
summary['First min Count umis per Cell']  = params.minCountPerCell1
summary['Second min Count umis per Cell']  = params.minCountPerCell2
summary['First min Count genes per Cell']  = params.minCountPerCellGene1
summary['Second min Count genes per Cell']  = params.minCountPerCellGene2
/*--------------------------------*/
summary['Annotation']   = params.genomeAnnotationPath
summary['Max Memory']     = params.maxMemory
summary['Max CPUs']       = params.maxCpus
summary['Max Time']       = params.maxTime
summary['Container Engine'] = workflow.containerEngine
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Config Profile'] = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


/*##########################   STEP 1: MAPPING  ####################################*/

process umiExtraction {
  tag "${prefix}"
  label 'umiTools'
  label 'highCpu'
  label 'highMem'
  publishDir "${params.outdir}/umiExtraction", mode: 'copy'

  input: 
  set val(prefix), file(reads) from rawReadsFastqc

  output:
  set val(prefix), file("*_UMIsExtracted.R1.fastq"), file("*_UMIsExtracted.R2.fastq") into chUmiExtracted
  set val(prefix), file("*_umiExtract.log") into chUmiExtractedLog
  file("v_umi_tools.txt") into chUmiToolsVersion

  script:
  length_umi = params.umi_size
  opts ="--extract-method=regex --bc-pattern='(?P<discard_1>.*ATTGCGCAATG)(?P<umi_1>.{$length_umi})(?P<discard_2>GGG).*' --stdin=${reads[0]} --stdout=${prefix}_UMIsExtracted.R1.fastq --read2-in=${reads[1]} --read2-out=${prefix}_UMIsExtracted.R2.fastq "
  """
  # Extract barcdoes and UMIs and add to read names
  umi_tools extract $opts --log=${prefix}_umiExtract.log 

  umi_tools --version &> v_umi_tools.txt
  """
}

process mergeReads {
  tag "${prefix}"
  label 'seqkit'
  label 'medCpu'
  label 'medMem'
  publishDir "${params.outdir}/mergeReads", mode: 'copy'

  input:
  set val(prefix), file(reads), file(umiReads_R1), file(umiReads_R2) from chMergeReadsFastq.join(chUmiExtracted)
  //set val(prefix), file(umiReads_R1), file(umiReads_R2) from chUmiExtracted

  output:
  set val(prefix), file("*_totReads.R1.fastq"), file("*_totReads.R2.fastq") into chMergeReads
  set val(prefix), file("*_umisReadsIDs.txt") into chUmiReadsIDs
  set val(prefix), file("*_NonUmisReadsIDs.txt") into chNonUmiReadsIDs
  set val(prefix), file("*_pUMIs.txt") into chCountSummaryExtUMI
  file("v_seqkit.txt") into chSeqkitVersion

  script:
  """
  # Get UMI read IDs (without UMIs)
  seqkit seq -n -i ${umiReads_R1}  > ${prefix}_umisReadsIDs.txt

  # Extract non umis reads
  seqkit grep -v -f ${prefix}_umisReadsIDs.txt ${reads[0]} -o ${prefix}_nonUMIs.R1.fastq
  seqkit grep -v -f ${prefix}_umisReadsIDs.txt ${reads[1]} -o ${prefix}_nonUMIs.R2.fastq

  # Get non UMI reads IDs
  seqkit seq -n -i ${prefix}_nonUMIs.R1.fastq > ${prefix}_NonUmisReadsIDs.txt

  # Merge non umis reads + umi reads (with umi in read names)
  cat ${umiReads_R1} > ${prefix}_totReads.R1.fastq
  cat ${prefix}_nonUMIs.R1.fastq >> ${prefix}_totReads.R1.fastq

  cat ${umiReads_R2} > ${prefix}_totReads.R2.fastq
  cat ${prefix}_nonUMIs.R2.fastq >> ${prefix}_totReads.R2.fastq

  ## Save % UMIs reads 
  nb_lines=`wc -l < <(gzip -cd ${reads[0]})`
  nb_totreads=\$(( \$nb_lines / 4 ))
  nb_umis=`wc -l < ${prefix}_umisReadsIDs.txt`
  echo "percentUMI:\$(( \$nb_umis * 100 / \$nb_totreads ))" > ${prefix}_pUMIs.txt

  # Get version 
  echo SeqKit > name
  seqkit --help | grep Version > version
  paste name version > v_seqkit.txt
  """
}

process trimReads{
  tag "${prefix}"
  label 'cutadapt'
  label 'medCpu'
  label 'medMem'
  publishDir "${params.outdir}/trimReads", mode: 'copy'

  input:
  set val(prefix), file(totReadsR1), file(totReadsR2) from chMergeReads

  output:
  set val(prefix), file("*_trimmed.R1.fastq"), file("*_trimmed.R2.fastq") into chTrimmedReads
  set val(prefix), file("*_trimmed.log") into chtrimmedReadsLog
  file("v_cutadapt.txt") into chCutadaptVersion

  script:
  """
  cutadapt -G XGCATACGAT{30} --minimum-length=15 -o ${prefix}_trimmed.R1.fastq -p ${prefix}_trimmed.R2.fastq ${totReadsR1} ${totReadsR2} > ${prefix}_trimmed.log
  cutadapt --version &> v_cutadapt.txt
  """
}


process readAlignment {
  tag "${prefix}"
  label 'STAR'
  label 'extraCpu'
  label 'extraMem'
  publishDir "${params.outdir}/readAlignment", mode: 'copy'

  input :
  file genomeIndex from chStar.collect()
  file genomeGtf from chGtfSTAR.collect()
  set val(prefix), file(trimmedR1) , file(trimmedR2) from chTrimmedReads
	
  output :
  set val(prefix), file("*Aligned.sortedByCoord.out.bam") into chAlignedBam
  file "*.out" into chAlignmentLogs
  file("v_star.txt") into chStarVersion

  script:  
  """
  STAR \
    --genomeDir $genomeIndex \
    --sjdbGTFfile $genomeGtf \
    --readFilesIn ${trimmedR1},${trimmedR2} \
    --runThreadN ${task.cpus} \
    --outFilterMultimapNmax 1 \
    --outFileNamePrefix ${prefix} \
    --outSAMtype BAM SortedByCoordinate 
    #--clip3pAdapterSeq CTGTCTCTTATACACATCT \
    #--limitSjdbInsertNsj 2000000 \
    #--outFilterIntronMotifs RemoveNoncanonicalUnannotated 

    # clip3pAdapterSeq = coupe l'adaptater 3' des R2 (~2%) 
    # limitSjdbInsertNsj = augmente le nombre de splice junctions à inserer
    # outFilterIntronMotifs = supprime les sp non annotées

  STAR --version &> v_star.txt
  """
}

process readsAssignment {
  tag "${prefix}"
  label 'featureCounts'
  label 'medCpu'
  label 'medMem'
  publishDir "${params.outdir}/readsAssignment", mode: 'copy'

  input :
  set val(prefix), file(alignedBam) from chAlignedBam
  file(genome) from chGtfFC.collect()

  output : 
  set val(prefix), file("*featureCounts.bam") into chAssignBam
  file "*.summary" into chAssignmentLogs
  file("v_featurecounts.txt") into chFCversion

  script:
  """	
  featureCounts  -p \
    -a ${genome} \
    -o ${prefix}_counts \
    -T ${task.cpus} \
    -R BAM \
    -g gene_name \
    ${alignedBam}

  featureCounts -v &> v_featurecounts.txt
  """
}

process sortBam {
  tag "${prefix}"
  label 'samtools'
  label 'medCpu'
  label 'medMem'
  publishDir "${params.outdir}/sortBam", mode: 'copy'

  input:
  set val(prefix), file(assignBam) from chAssignBam
	
  output:
  set val(prefix), file("*_Sorted.bam") into chSortedBAM_bigWig, chSortedBAM_sepReads, chSortedBAM_readCounts
  file("v_samtools.txt") into chSamtoolsVersion

  script :
  """
  samtools sort -@ ${task.cpus} ${assignBam} -o ${prefix}_Sorted.bam

  samtools --version &> v_samtools.txt
  """
}

process separateReads {
  tag "${prefix}"
  label 'samtools'
  label 'medCpu'
  label 'medMem'
  publishDir "${params.outdir}/separateReads", mode: 'copy'

  input :
  set val(prefix), file(sortedBam), file(umisReadsIDs), file(nonUmisReadsIDs) from chSortedBAM_sepReads.join(chUmiReadsIDs).join(chNonUmiReadsIDs)
  //set val(prefix), file(umisReadsIDs) from chUmiReadsIDs
  //set val(prefix), file(nonUmisReadsIDs) from chNonUmiReadsIDs

  output:
  set val("${prefix}_umi"), file("*_assignedUMIs.bam") into chUmiBam, chUmiBam_countMtx
  set val("${prefix}_NonUmi"), file("*_assignedNonUMIs.bam") into chNonUmiBam

  script:  
  """
  # Separate umi and non umi reads
  samtools view ${sortedBam} > ${prefix}assignedAll.sam

  # save header and extract umi reads 
  samtools view -H ${sortedBam} > ${prefix}_assignedUMIs.sam
  fgrep -f ${umisReadsIDs} ${prefix}assignedAll.sam >> ${prefix}_assignedUMIs.sam
  # sam to bam
  samtools view -bh ${prefix}_assignedUMIs.sam > ${prefix}_assignedUMIs.bam

  # save header and extract non umi reads 
  samtools view -H ${sortedBam} > ${prefix}_assignedNonUMIs.sam
  fgrep -f ${nonUmisReadsIDs} ${prefix}assignedAll.sam >> ${prefix}_assignedNonUMIs.sam
  # sam to bam
  samtools view -bh ${prefix}_assignedNonUMIs.sam > ${prefix}_assignedNonUMIs.bam

  #-h Include the header in the output.
  #-H Output the header only.
  """
}

process countMatrices {
  tag "${prefix}"
  label 'umiTools'
  label 'medCpu'
  label 'medMem'
  publishDir "${params.outdir}/countMatrices", mode: 'copy'

  input:
  set val(prefix), file(umiBam) from chUmiBam_countMtx

  output:
  set val(prefix), file("*_Counts.tsv.gz") into chMatrices, chMatrices_dist, chMatrices_counts
  set val(prefix), file("*_UmiCounts.log") into chMatricesLog

  script:
  """
  # Count UMIs per gene per cell
  samtools index ${umiBam}
  umi_tools count --method=cluster --per-gene --gene-tag=XT --assigned-status-tag=XS -I ${umiBam} -S ${prefix}_Counts.tsv.gz > ${prefix}_UmiCounts.log
  """
}

process bigWig {
  tag "${prefix}"
  label 'deeptools'
  label 'extraCpu'
  label 'extraMem'
  publishDir "${params.outdir}/bigWig", mode: 'copy'

  input:
  //set val(prefix), file(assignedBam) from chSortedBAM_bigWig
  set val(prefix), file(bam) from chSortedBAM_bigWig.concat(chUmiBam).concat(chNonUmiBam) //  _NonUmi_assignedNonUMIs.bam _umi_assignedNonUMIs.bam

  output:
  set val(prefix), file("*_coverage.bw") into chBigWig // L386_coverage.bw , L386_umi_coverage.bw, L386_NonUmi_coverage.bw
  set val(prefix), file("*_coverage.log") into chBigWigLog
  file("v_deeptools.txt") into chBamCoverageVersion

  script:
  """
  ## Create bigWig files
  samtools index ${bam}
  #bamCoverage --normalizeUsing CPM -b ${bam} -of bigwig -o ${prefix}_coverage.bw > ${prefix}_coverage.log
  bamCoverage -b ${bam} -of bigwig -o ${prefix}_coverage.bw --numberOfProcessors=5 > ${prefix}_coverage.log

  bamCoverage --version &> v_deeptools.txt
  """
}

//chUmiBam.view()

process genebody_coverage {
  tag "${prefix}"
  label 'rseqc'
  label 'extraCpu'
  label 'extraMem'
  
  publishDir "${params.outdir}/genebody_coverage" , mode: 'copy',
  saveAs: {filename ->
      if (filename.indexOf("geneBodyCoverage.curves.pdf") > 0)       "geneBodyCoverage/$filename"
      else if (filename.indexOf("geneBodyCoverage.r") > 0)           "geneBodyCoverage/rscripts/$filename"
      else if (filename.indexOf("geneBodyCoverage.txt") > 0)         "geneBodyCoverage/data/$filename"
      else if (filename.indexOf("log.txt") > -1) false
      else filename
  }

  input:
  //set val(prefix), file (bigwig) from chBigWig.filter( ~/.*mi_coverage.bw/ ) // L386_coverage.bw, L386_umi_coverage.bw, L386_NonUmi_coverage.bw
  set val(prefix), file (bg) from chBigWig.filter( ~/.*mi_*/ )
  //set val(prefix), file (bg) from chBigWig.filter( ~/.*mi_.*/ )
  file bed12 from chBedGeneCov.collect()
  // channel = pile
  // quand site tous les fichiers => c'est que commandes differentes sur les deux
  // To create new channel from one, those == spikego in 1 
  // chAlignReads.choice( chAlignSpike, chAlignRef ){ it -> it[1] =~ 'spike' ? 1 : 0 }

  output:
  file "*.{txt,pdf,r}" into chGeneCov_res

  script:
  """
  #samtools index ${bg}
  geneBody_coverage2.py \\
      -i ${bg} \\
      -o ${prefix}.rseqc \\
      -r $bed12

  #geneBody_coverage.py \\
  #    -i ${bam} \\
  #    -o ${prefix}.rseqc \\
  #    -r $bed12
  #mv log.txt ${prefix}.rseqc.log.txt
  """
}

//chBigWig.view()


/*##########################   STEP 2: CELL VIABILITY  ####################################*/

process umiPerGeneDist{
  tag "${prefix}"
  label 'R'
  label 'lowCpu'
  label 'lowMem'
  publishDir "${params.outdir}/umiPerGeneDist", mode: 'copy'

  input:
  set val(prefix), file(matrix) from chMatrices_dist

  output:
  set val(prefix), file ("*_HistUMIperGene_mqc.csv") into chUMIperGene

  script:
  """
  umiPerGene_dist.r ${matrix} ${prefix}
  """ 
}

process countUMIGenePerCell{
  tag "${prefix}"
  label 'R'
  label 'lowCpu'
  label 'lowMem'
  publishDir "${params.outdir}/countUMIGenePerCell", mode: 'copy'

  input:
  set val(prefix), file(matrix) from chMatrices_counts

  output:
  //set val(prefix), file ("*_nbGenePerCell_mqc.csv") into chGenePerCell
  //set val(prefix), file ("*_nbUMIPerCell_mqc.csv") into chUmiPerCell
  set val(prefix), file ("*_countPerCell_mqc.csv") into chUMI_Gene_perCell

  script:
  """
  umiGenePerCell.r ${matrix} ${prefix}
  """ 
}

process cellAnalysis{
  tag "${prefix}"
  label 'R'
  label 'highCpu'
  label 'highMem'
  publishDir "${params.outdir}/cellAnalysis", mode: 'copy'

  input:
  file ("matrices/${prefix}*") from chMatrices.collect()

  output:
  file ("10Xoutput/") into ch10X
  file ("resume.txt") into chResume
  //file ("HistUMIperGene.mqc") into chUMIperGene
  //file("UMIGenesPerCell_mqc.csv") into chUMI_Gene_perCell
  //file ("HistUMIperCell_mqc.csv") into chUMIperCell
  //file ("HistGenePerCell_mqc.csv") into chGenesPerCell
  file ("RatioPerCell_mqc.csv") into chUmiGeneRatio
  file ("MtGenePerCell_mqc.csv") into chMT
  file ("v_R.txt") into chRversion

  script:
  """
  cellViability.r matrices/ 10Xoutput/
  R --version &> v_R.txt  
  """ 
}

/*-----------------------------------*/

/***********
 * MultiQC *
 ***********/

process getSoftwareVersions{
  label 'python'
  label 'lowCpu'
  label 'lowMem'
  publishDir path: "${params.outdir}/software_versions", mode: "copy"

  when:
  !params.skipSoftVersions

  input:
  file("v_umi_tools.txt") from chUmiToolsVersion.first().ifEmpty([])
  file("v_seqkit.txt") from chSeqkitVersion.first().ifEmpty([])
  file("v_cutadapt.txt") from chCutadaptVersion.first().ifEmpty([])
  file("v_star.txt") from chStarVersion.first().ifEmpty([])
  file("v_featurecounts.txt") from chFCversion.first().ifEmpty([])
  file("v_samtools.txt") from chSamtoolsVersion.first().ifEmpty([])
  file("v_deeptools.txt") from chBamCoverageVersion.first().ifEmpty([])
  file ("v_R.txt") from chRversion.ifEmpty([])

  output:
  file 'software_versions_mqc.yaml' into softwareVersionsYaml

  script:
  """
  echo $workflow.manifest.version &> v_pipeline.txt
  echo $workflow.nextflow.version &> v_nextflow.txt
  scrape_software_versions.py &> software_versions_mqc.yaml
  """
}

process workflowSummaryMqc {
  when:
  !params.skipMultiQC

  output:
  file 'workflow_summary_mqc.yaml' into workflowSummaryYaml

  exec:
  def yaml_file = task.workDir.resolve('workflow_summary_mqc.yaml')
  yaml_file.text  = """
  id: 'summary'
  description: " - this information is collected when the pipeline is started."
  section_name: 'Workflow Summary'
  section_href: 'https://gitlab.curie.fr/data-analysis/chip-seq'
  plot_type: 'html'
  data: |
        <dl class=\"dl-horizontal\">
  ${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
  """.stripIndent()
}

process multiqc {
  label 'multiqc'
  label 'medCpu'
  label 'medMem'
  publishDir "${params.outdir}/MultiQC", mode: 'copy'

  when:
  !params.skipMultiQC

  input:
  file splan from chSplan.collect()
  file multiqcConfig from chMultiqcConfig
  file design from chDesignMqc.collect().ifEmpty([])
  file metadata from chMetadata.ifEmpty([])
  file ('software_versions/*') from softwareVersionsYaml.collect().ifEmpty([])
  file ('workflow_summary/*') from workflowSummaryYaml.collect()
  //Modules
  file ('trimming/*') from chtrimmedReadsLog.collect()
  file ('star/*') from chAlignmentLogs.collect()
  file ('FC/*') from chAssignmentLogs.collect()
  file ('coverage/*') from chGeneCov_res.collect().ifEmpty([])
  //LOGS
  file ('umiExtract/*') from chUmiExtractedLog.collect()
  file('mergeReads/*') from chCountSummaryExtUMI.collect()
  file ('bigwig/*') from chBigWigLog.collect()
  file (resume) from chResume
  //PLOTS
  file ("umiPerGene/*") from chUMIperGene.collect() // linegraph == histogram
  //file ("nbUMI/*") from chUmiPerCell.collect()  // bargraph
  //file ("nbGene/*") from chGenePerCell.collect() // bargraph 
  file (umi_Gene_perCell) from chUMI_Gene_perCell //
  file ("ratio/*") from chUmiGeneRatio.collect() // UmiGenePerCell_mqc.csv
  file ("mt/*") from chMT.collect() // MtGenePerCell_mqc.csv

  output: 
  file splan
  file "*report.html" into multiqc_report
  file "*_data"

  script:
  rtitle = customRunName ? "--title \"$customRunName\"" : ''
  rfilename = customRunName ? "--filename " + customRunName + "_report" : "--filename report"
  metadataOpts = params.metadata ? "--metadata ${metadata}" : ""
  //isPE = params.singleEnd ? "" : "-p"
  designOpts= params.design ? "-d ${params.design}" : ""
  modules_list = "-m custom_content -m cutadapt -m samtools -m star -m featureCounts -m deeptools  -m rseqc"

  //umisFiltre1 = params.minCountPerCell1 ? "--minCountPerCell1 ${params.minCountPerCell1}" : ""
  //umisFiltre2 = params.minCountPerCell2 ? "--minCountPerCell2 ${params.minCountPerCell2}" : ""
  //umisFiltreGene1 = params.minCountPerCellGene1 ? "--minCountPerCellGene1 ${params.minCountPerCellGene1}" : ""
  //umisFiltreGene2 = params.minCountPerCellGene2 ? "--minCountPerCellGene2 ${params.minCountPerCellGene2}" : ""
  // ${umisFiltre1} ${umisFiltre2} ${umisFiltreGene1} ${umisFiltreGene2}

  """
  stat2mqc.sh ${splan} 
  mqc_header.py --splan ${splan} --name "PIPELINE" --version ${workflow.manifest.version} ${metadataOpts} > multiqc-config-header.yaml
  multiqc . -f $rtitle $rfilename -c multiqc-config-header.yaml -c $multiqcConfig $modules_list
  """
}

/****************
 * Sub-routines *
 ****************/

process outputDocumentation {
  label 'python'
  label 'lowCpu'
  label 'lowMem'

  publishDir "${params.outdir}/pipeline_info", mode: 'copy'

  input:
  file output_docs from chOutputDocs
  file images from chOutputDocsImages

  output:
  file "results_description.html"

  script:
  """
  markdown_to_html.py $output_docs -o results_description.html
  """
}

workflow.onComplete {

  // pipelineReport.html
  def reportFields = [:]
  reportFields['version'] = workflow.manifest.version
  reportFields['runName'] = customRunName ?: workflow.runName
  reportFields['success'] = workflow.success
  reportFields['dateComplete'] = workflow.complete
  reportFields['duration'] = workflow.duration
  reportFields['exitStatus'] = workflow.exitStatus
  reportFields['errorMessage'] = (workflow.errorMessage ?: 'None')
  reportFields['errorReport'] = (workflow.errorReport ?: 'None')
  reportFields['commandLine'] = workflow.commandLine
  reportFields['projectDir'] = workflow.projectDir
  reportFields['summary'] = summary
  reportFields['summary']['Date Started'] = workflow.start
  reportFields['summary']['Date Completed'] = workflow.complete
  reportFields['summary']['Pipeline script file path'] = workflow.scriptFile
  reportFields['summary']['Pipeline script hash ID'] = workflow.scriptId
  if(workflow.repository) reportFields['summary']['Pipeline repository Git URL'] = workflow.repository
  if(workflow.commitId) reportFields['summary']['Pipeline repository Git Commit'] = workflow.commitId
  if(workflow.revision) reportFields['summary']['Pipeline Git branch/tag'] = workflow.revision

  // Render the TXT template
  def engine = new groovy.text.GStringTemplateEngine()
  def tf = new File("$baseDir/assets/onCompleteTemplate.txt")
  def txtTemplate = engine.createTemplate(tf).make(reportFields)
  def reportTxt = txtTemplate.toString()

  // Render the HTML template
  def hf = new File("$baseDir/assets/onCompleteTemplate.html")
  def htmlTemplate = engine.createTemplate(hf).make(reportFields)
  def reportHtml = htmlTemplate.toString()

  // Write summary HTML to a file
  def outputSummaryDir = new File( "${params.summaryDir}/" )
  if( !outputSummaryDir.exists() ) {
    outputSummaryDir.mkdirs()
  }
  def outputHtmlFile = new File( outputSummaryDir, "pipelineReport.html" )
  outputHtmlFile.withWriter { w -> w << reportHtml }
  def outputTxtFile = new File( outputSummaryDir, "pipelineReport.txt" )
  outputTxtFile.withWriter { w -> w << reportTxt }

  // onComplete file
  File woc = new File("${params.outdir}/onComplete.txt")
  Map endSummary = [:]
  endSummary['Completed on'] = workflow.complete
  endSummary['Duration']     = workflow.duration
  endSummary['Success']      = workflow.success
  endSummary['exit status']  = workflow.exitStatus
  endSummary['Error report'] = workflow.errorReport ?: '-'
  String endWfSummary = endSummary.collect { k,v -> "${k.padRight(30, '.')}: $v" }.join("\n")
  println endWfSummary
  String execInfo = "Execution summary\n${endWfSummary}\n"
  woc.write(execInfo)

  // final logs
  if(workflow.success){
      log.info "Pipeline Complete"
  }else{
    log.info "FAILED: $workflow.runName"
  }
}
