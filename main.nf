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
devMessageFile = file("$projectDir/assets/devMessage.txt")

def helpMessage() {
  if ("${workflow.manifest.version}" =~ /dev/ ){
    devMess = file("$projectDir/assets/devMessage.txt")
    log.info devMessageFile.text
  }

  log.info """
  SmartSeq3 v${workflow.manifest.version}
  ======================================================================

  Usage:
  nextflow run main.nf --reads '*_R{1,2}.fastq.gz' --genome 'hg38' -profile conda
  nextflow run main.nf --samplePlan samplePlan --genome 'hg38' -profile conda

  Mandatory arguments:
    --reads [file]                Path to input data (must be surrounded with quotes)
    --samplePlan [file]           Path to sample plan input file (cannot be used with --reads)
    --genome [str]                Name of genome reference
    -profile [str]                Configuration profile to use. test / conda / multiconda / path / multipath / singularity / docker / cluster (see below)
    --protocol [str]
  
  Inputs:
    --starIndex [dir]             Index for STAR aligner
    --gtf [dir]                   GTF for STAR aligner

  Skip options: All are false by default
    --skipSoftVersion [bool]      Do not report software version. Default is false.
    --skipMultiQC [bool]          Skips MultiQC. Default is false.
    --skipGeneCov [bool]          Skips calculating genebody coverage. Default is false.
    --skipSatCurves [bool]        Skips saturation curves. Default is false.
  
  Genomes: If not specified in the configuration file or if you wish to overwrite any of the references given by the --genome field
  --genomeAnnotationPath [file]      Path  to genome annotation folder

  Other options:
    --outDir [file]               The output directory where the results will be saved
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
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
  exit 1, "The provided genome '${params.genome}' is not available in the genomes.config file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

//Stage config files
Channel
  .fromPath(params.multiqcConfig, checkIfExists: true)
  .set{chMultiqcConfig}
chOutputDocs = file("$projectDir/docs/output.md", checkIfExists: true)
chOutputDocsImages = file("$projectDir/docs/images/", checkIfExists: true)

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
    .into { chGtfSTAR; chGtfFCumis; chGtfFCallreads }
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

// Create a channel for input read files
if(params.samplePlan){
  if(params.singleEnd){
    Channel
      .from(file("${params.samplePlan}"))
      .splitCsv(header: false)
      .map{ row -> [ row[0], [file(row[2])]] }
      .into { rawReads }
  }else{
    Channel
      .from(file("${params.samplePlan}"))
      .splitCsv(header: false)
      .map{ row -> [ row[0], [file(row[2]), file(row[3])]] }
      .into { rawReads }
   }
  params.reads=false
}
else if(params.readPaths){
  if(params.singleEnd){
    Channel
      .from(params.readPaths)
      .map { row -> [ row[0], [file(row[1][0])]] }
      .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied." }
      .into { rawReads }
  } else {
    Channel
      .from(params.readPaths)
      .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
      .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied." }
      .into { rawReads }
  }
} else {
  Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .into { rawReads }
}

// Make sample plan if not available
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
summary['Pipeline Name']  = 'SmartSeq3'
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = customRunName ?: workflow.runName
summary['Command Line'] = workflow.commandLine
if (params.samplePlan) {
   summary['SamplePlan']   = params.samplePlan
}else{
   summary['Reads']        = params.reads
}
summary['Genome']       = params.genome
summary['Annotation']   = params.genomeAnnotationPath
summary['Max Memory']     = params.maxMemory
summary['Max CPUs']       = params.maxCpus
summary['Max Time']       = params.maxTime
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outDir
summary['Config Profile'] = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

process umiExtraction {
  tag "${prefix}"
  label 'umiTools'
  label 'medCpu'
  label 'medMem'

  publishDir "${params.outDir}/umiExtraction", mode: 'copy'

  input: 
  set val(prefix), file(reads) from rawReads

  output:
  set val(prefix), file("*_totReads.R1.fastq.gz"), file("*_totReads.R2.fastq.gz") into chTrimTag
  set val(prefix), file("*_umiExtractR1.log"), file("*_umiExtractR2.log") into chUmiExtractedLog
  set val(prefix), file("*_nonUmisReadsIDs.txt") into chUmiReadsIDs_exUMIreads, chUmiReadsIDs_exNonUMIreads
  set val(prefix), file("*_pUMIs.txt") into chCountSummaryExtUMI
  set val(prefix), file("*_nbTotFrag.txt") into chTotFrag

  file("v_umi_tools.txt") into chUmiToolsVersion

  script:
  """
  # Extract sequences that have tag+UMI+GGG and add UMI to read names
  # If no umi is find, the reads is leave without changement

  if [[ ${params.protocol} == "flashseq" ]]
  then 
    #tag="AAGCAGTGGTATCAACGCAGAGT"
    umi_tools extract --extract-method=regex --bc-pattern='(?P<discard_1>.*CGCAGAGT)(?P<umi_1>.{$params.umi_size})(?P<discard_2>GGG).*' \\
                    --stdin=${reads[0]} --stdout=${prefix}_UMIsExtractedR1.R1.fastq.gz \\
                    --read2-in=${reads[1]} --read2-out=${prefix}_UMIsExtractedR1.R2.fastq.gz \\
                    --filtered-out ${prefix}_noUMIinR1.R1.fastq.gz --filtered-out2 ${prefix}_noUMIinR1.R2.fastq.gz \\
                    --log=${prefix}_umiExtractR1.log

    # use R2 as stdin to exrtract umis also from R2 
    umi_tools extract --extract-method=regex --bc-pattern='(?P<discard_1>.*CGCAGAGT)(?P<umi_1>.{$params.umi_size})(?P<discard_2>GGG).*' \\
                    --stdin=${prefix}_noUMIinR1.R2.fastq.gz  --stdout=${prefix}_UMIsExtractedR2.R2.fastq.gz \\
                    --read2-in=${prefix}_noUMIinR1.R1.fastq.gz  --read2-out=${prefix}_UMIsExtractedR2.R1.fastq.gz \\
                    --filtered-out ${prefix}_nonUMIreads.R2.fastq.gz --filtered-out2 ${prefix}_nonUMIreads.R1.fastq.gz \\
                    --log=${prefix}_umiExtractR2.log                
  else 
    #tag="TGCGCAATG"
    # extract umis from R1 
    umi_tools extract --extract-method=regex --bc-pattern='(?P<discard_1>.*TGCGCAATG)(?P<umi_1>.{$params.umi_size})(?P<discard_2>GGG).*' \\
                    --stdin=${reads[0]} --stdout=${prefix}_UMIsExtractedR1.R1.fastq.gz \\
                    --read2-in=${reads[1]} --read2-out=${prefix}_UMIsExtractedR1.R2.fastq.gz \\
                    --filtered-out ${prefix}_noUMIinR1.R1.fastq.gz --filtered-out2 ${prefix}_noUMIinR1.R2.fastq.gz \\
                    --log=${prefix}_umiExtractR1.log
    # extract umi from R2 
    umi_tools extract --extract-method=regex --bc-pattern='(?P<discard_1>.*TGCGCAATG)(?P<umi_1>.{$params.umi_size})(?P<discard_2>GGG).*' \\
                    --stdin=${prefix}_noUMIinR1.R2.fastq.gz  --stdout=${prefix}_UMIsExtractedR2.R2.fastq.gz \\
                    --read2-in=${prefix}_noUMIinR1.R1.fastq.gz  --read2-out=${prefix}_UMIsExtractedR2.R1.fastq.gz \\
                    --filtered-out ${prefix}_nonUMIreads.R2.fastq.gz --filtered-out2 ${prefix}_nonUMIreads.R1.fastq.gz \\
                    --log=${prefix}_umiExtractR2.log   
  fi

  # save read IDs that have no umis in a file to extract them after alignment 
  seqkit seq -j ${task.cpus} -n -i ${prefix}_nonUMIreads.R2.fastq.gz -o ${prefix}_nonUmisReadsIDs.txt

  # concatenate R1 and R2 umi reads == all umi reads 
  cat ${prefix}_UMIsExtractedR2.R1.fastq.gz >> ${prefix}_UMIsExtractedR1.R1.fastq.gz 
  ############## Save % UMIs reads
  nb_lines=`wc -l < <(gzip -cd ${reads[0]})`
  nb_totFrag=\$(( \$nb_lines / 4 ))
  echo "totFrag: \$nb_totFrag" > ${prefix}_nbTotFrag.txt

  nb_lines=`wc -l < <(gzip -cd ${prefix}_UMIsExtractedR1.R1.fastq.gz) `
  nb_umis=\$(( \$nb_lines / 4 ))
  echo "percentUMI:\$(( \$nb_umis * 100 / \$nb_totFrag ))" > ${prefix}_pUMIs.txt
  ##############
  # add non umi reads == all reads 
  cat ${prefix}_nonUMIreads.R1.fastq.gz >> ${prefix}_UMIsExtractedR1.R1.fastq.gz 
  mv ${prefix}_UMIsExtractedR1.R1.fastq.gz ${prefix}_totReads.R1.fastq.gz

  # concatenate R1 and R2 umi reads
  cat ${prefix}_UMIsExtractedR2.R2.fastq.gz >> ${prefix}_UMIsExtractedR1.R2.fastq.gz 
  # add non umi reads
  cat ${prefix}_nonUMIreads.R2.fastq.gz >> ${prefix}_UMIsExtractedR1.R2.fastq.gz 
  mv ${prefix}_UMIsExtractedR1.R2.fastq.gz ${prefix}_totReads.R2.fastq.gz

  umi_tools --version &> v_umi_tools.txt
  """
}

process trimTag{
  tag "${prefix}"
  label 'cutadapt'
  label 'highCpu'
  label 'highMem'

  publishDir "${params.outDir}/trimTag", mode: 'copy'

  input:
  set val(prefix), file(totReadsR1), file(totReadsR2) from chTrimTag

  output:
  set val(prefix), file("*_trimmed.R1.fastq.gz"), file("*_trimmed.R2.fastq.gz") into chTrimmedReads 
  set val(prefix), file("*_trimSens.log"), file("*_trimAntisens.log") into chtrimmedReadsLog
  file("v_cutadapt.txt") into chCutadaptVersion

  script:
  """
  ## Delete linker + polyA/T queue
  ## cutadpat look for the tag+polyA/T in 5' and cut them if some bases are found 

  if [[ ${params.protocol} == "flashseq" ]]
  then 
    # 1) sens strand = AAAAAtag
    cutadapt \
    -g A{30}GTACTCTGCGTTGATACCACTGCTT -g A{20}GTACTCTGCGTTGATACCACTGCTT -g A{15}GTACTCTGCGTTGATACCACTGCTT \
    -G A{30}GTACTCTGCGTTGATACCACTGCTT -G A{20}GTACTCTGCGTTGATACCACTGCTT -G A{15}GTACTCTGCGTTGATACCACTGCTT \
    --minimum-length=20 \
    --cores=${task.cpus} -o ${prefix}_trimSens.R1.fastq.gz -p ${prefix}_trimSens.R2.fastq.gz \
    <(gzip -cd ${totReadsR1}) <(gzip -cd ${totReadsR2}) > ${prefix}_trimSens.log

    # 2) antisens strand == check in 5' part of my R1 (-g) and my R2 (-G) 
  
    cutadapt \
    -g AAGCAGTGGTATCAACGCAGAGTACT{30} -g AAGCAGTGGTATCAACGCAGAGTACT{20} -g AAGCAGTGGTATCAACGCAGAGTACT{15} \
    -G AAGCAGTGGTATCAACGCAGAGTACT{30} -G AAGCAGTGGTATCAACGCAGAGTACT{20} -G AAGCAGTGGTATCAACGCAGAGTACT{15} \
    --minimum-length=20 \
    --cores=${task.cpus} \
    -o ${prefix}_trimmed.R1.fastq.gz -p ${prefix}_trimmed.R2.fastq.gz \
    ${prefix}_trimSens.R1.fastq.gz ${prefix}_trimSens.R2.fastq.gz > ${prefix}_trimAntisens.log

  else 

    #5'- cDNA(pA)TCGTATGCTGCT GATGCTCGT -3' 
    #3'- cDNA(dT)AGCATACGACGA CTACGAGCA -5' 

    # 1) sens strand
    # -a remove matches between the polyA and the 3' end
    cutadapt -a A{30}TCGTATGCTGCT -a A{20}TCGTATGCTGCT \
    -A A{30}TCGTATGCTGCT -A A{20}TCGTATGCTGCT \
    --minimum-length=20 \
    --cores=${task.cpus} -o ${prefix}_trimSens.R1.fastq.gz -p ${prefix}_trimSens.R2.fastq.gz \
    ${totReadsR1} ${totReadsR2} > ${prefix}_trimSens.log
    
    # 2) antisens strand 
    # echo AGCATACGACGACTACGAGCA | rev = ACGAGCATCAGCAGCATACGA
    cutadapt -g ATCAGCAGCATACGAT{30} -g ATCAGCAGCATACGAT{20} \
    -G ATCAGCAGCATACGAT{30} -G ATCAGCAGCATACGAT{20} \
    --minimum-length=20 \
    --cores=${task.cpus} -o ${prefix}_trimmed.R1.fastq.gz -p ${prefix}_trimmed.R2.fastq.gz \
    ${prefix}_trimSens.R1.fastq.gz ${prefix}_trimSens.R2.fastq.gz > ${prefix}_trimAntisens.log

  fi
  
  cutadapt --version &> v_cutadapt.txt
  """
}

// From nf-core
// Function that checks the alignment rate of the STAR output
// and returns true if the alignment passed and otherwise false
skippedPoorAlignment = []
def checkStarLog(logs) {
  def percentAligned = 0;
  logs.eachLine { line ->
    if ((matcher = line =~ /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/)) {
      percentAligned = matcher[0][1]
    }else if ((matcher = line =~ /Uniquely mapped reads number\s*\|\s*([\d\.]+)/)) {
      numAligned = matcher[0][1]
    }
  }
  logname = logs.getBaseName() - 'Log.final'
  if(numAligned.toInteger() <= 2000.toInteger() ){
      log.info "#################### LESS THAN 2000 READS! IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! ($logname)  >> ${percentAligned}% <<"
      skippedPoorAlignment << logname
      return false
  } else {
      log.info "          Passed alignment > star ($logname)   >> ${percentAligned}% <<"
      return true
  }
}

// Update input channel
chStarRawReads = Channel.empty()
chStarRawReads = chTrimmedReads

process readAlignment {
  tag "${prefix}"
  label 'star'
  label 'extraCpu'
  label 'highMem'

  publishDir "${params.outDir}/star", mode: 'copy'

  input :
  file genomeIndex from chStar.collect()
  file genomeGtf from chGtfSTAR.collect()
  set val(prefix), file(trimmedR1) , file(trimmedR2) from chStarRawReads
	
  output :
  set val(prefix), file("*Log.final.out"), file("*Aligned.sortedByCoord.out.bam") into chAlignBam
  file "*.out" into chAlignmentLogs
  file("v_star.txt") into chStarVersion

  script:  
  """  
  STAR \
    --genomeDir $genomeIndex \
    --readFilesIn ${trimmedR1} ${trimmedR2} \
    --readFilesCommand zcat \
    --runThreadN ${task.cpus} \
    --outFilterMultimapNmax 1 \
    --outFileNamePrefix ${prefix} \
    --outSAMtype BAM SortedByCoordinate \
    --limitSjdbInsertNsj 2000000 \
    --sjdbGTFfile $genomeGtf --outFilterIntronMotifs RemoveNoncanonicalUnannotated 

    # outFilterMultimapNmax = max nb of loci the read is allowed to map to. If more, the read is concidered "map to too many loci". 
    # limitSjdbInsertNsj = max number of junctions to be insterted to the genome (those known (annotated) + those not annot. but found in many reads). 
    # Default is 1 000 000. By increasing it, more new junctions can be discovered. 
    # outFilterIntronMotifs = delete unannotated (not in genomeGtf) splice junctions with non-canonical (<=> unusual) intron motifs.
    # Non-canonical but annot. or canonical but not annot. will be kept.
    # NB: Canonical <=> juctions describe as having GT/AG, GC/AG or AT/AC (donor/acceptor) dinucleotide combination. 
    # Non-canonical are all other dinucleotide combinations. 

  STAR --version &> v_star.txt

  rm -rf ${prefix}_STARtmp ${prefix}_STARgenome
  """
}

// Filter removes all 'aligned' channels that fail the check
chAlignBam
  .filter { prefix, logs, bams -> checkStarLog(logs) }
  .map() {item -> [item[0], item[2]] }
  .dump (tag:'starbams')
  .into {chAlignBamCheck_umi; chAlignBamCheck_allreads}



// Umi Reads  ----------------------------------------------------------------------------------//

process extractUMIreads {
  tag "${prefix}"
  label 'samtools'
  label 'medCpu'
  label 'extraMem'

  publishDir "${params.outDir}/extractUMIreads", mode: 'copy'

  input :
  set val(prefix), file(aligned), file(nonUmisReadsIDs) from chAlignBamCheck_umi.join(chUmiReadsIDs_exUMIreads)

  output:
  set val("${prefix}"), file("*_Umi_assignedUMIs.{bam,bam.bai}") into chUmiBam, chUmiBam_geneBody
  set val(prefix), file("*_list_id_umi_seq.txt") into chListIDumis

  script:  
  """
  # Separate umi and non umi reads
  samtools view ${aligned} > ${prefix}assignedAll.sam

  # save header for extracting umi reads 
  samtools view -H ${aligned} > ${prefix}_assignedUMIs.sam

  nbLines=\$(wc -l < ${nonUmisReadsIDs})
  if((\$nbLines!=0))
  then
    fgrep -v -f ${nonUmisReadsIDs} ${prefix}assignedAll.sam >> ${prefix}_assignedUMIs.sam || echo "no sequence found"
  else
    cat ${prefix}assignedAll.sam >> ${prefix}_assignedUMIs.sam
  fi

  cut -f1 ${prefix}_assignedUMIs.sam > ${prefix}_list_id_umi_seq.txt

  # sam to bam
  samtools view -bh ${prefix}_assignedUMIs.sam > ${prefix}_Umi_assignedUMIs.bam

  # index
  samtools index ${prefix}_Umi_assignedUMIs.bam
  rm *.sam
  """
}

process readAssignmentUmis {
  tag "${prefix}" //${prefix}_Umi
  label 'featureCounts'
  label 'highCpu'
  label 'medMem'

  publishDir "${params.outDir}/featurecounts/umi", mode: 'copy'

  input :
  set val(prefix), file(aligned) from chUmiBam
  file(genome) from chGtfFCumis.collect()

  output : 
  set val(prefix), file("*featureCounts.bam") into chAssignBam
  file "*.summary" into chAssignmentLogs_umis
  set val(prefix), file("*_counts") into featureCountMatrix_umi
  file("v_featurecounts.txt") into chFCversion

  script:
  """	
  featureCounts  -p \
    -a ${genome} \
    -o ${prefix}_counts \
    -T ${task.cpus} \
    -R BAM \
    -g gene_name ${aligned[0]}

  featureCounts -v &> v_featurecounts.txt

  # -a annotation file
  # -R results format
  """
}

process sortAndIndexBam {
  tag "${prefix}" //${prefix}_Umi
  label 'samtools'
  label 'highCpu'
  label 'medMem'

  publishDir "${params.outDir}/bams/umi", mode: 'copy'

  input:
  set val(prefix), file(assignBam) from chAssignBam
	
  output:
  set val(prefix), file("*_Sorted.{bam,bam.bai}") into chSortedBAM, chAssignBam_dedup
  file("v_samtools.txt") into chSamtoolsVersion

  script :
  """
  samtools sort -@ ${task.cpus} ${assignBam} -o ${prefix}_Sorted.bam
  samtools index ${prefix}_Sorted.bam
  samtools --version &> v_samtools.txt
  """
}
// dedup and count in the same time
process countMatricesUMIs {
  tag "${prefix}" //${prefix}_Umi
  label 'umiTools'
  label 'medCpu'
  label 'medMem'

  publishDir "${params.outDir}/countMatricesUMI_logs", mode: 'copy'

  input:
  set val(prefix), file(umiBam) from chSortedBAM

  output:
  set val(prefix), file("*_Counts.tsv.gz") into chMatricesUMI, chMatrices_dist, chMatrices_UMIGenePerCell
  set val(prefix), file("*_UmiCounts.log") into chMatricesLog

  script:
  """
  # Count UMIs per gene per cell
  umi_tools count --method=cluster --per-gene --gene-tag=XT --assigned-status-tag=XS -I ${umiBam[0]} -S ${prefix}_Counts.tsv.gz > ${prefix}_UmiCounts.log
  """
}

process mergeUMIMatrices {
  label 'R'
  label 'highCpu'
  label 'medMem'

  publishDir "${params.outDir}/Matrices", mode: 'copy'

  input:
  file (umimatrices) from chMatricesUMI.collect()

  output:
  file ("10XlikeMatrix_umi/") into ch10X_mt
  file ("10XlikeMatrix_umi.zip") into ch10Xzip
  file ("UMI_gene_per_cell.txt") into chUmiResume, chUmiResume_mt
  file ("v_R.txt") into chRversion

  script:
  """
  merge_matrices.r 10XlikeMatrix_umi/ umi
  zip 10XlikeMatrix_umi.zip 10XlikeMatrix_umi/*
  R --version &> v_R.txt  
  """ 
}

// Non Umi reads  -------------------------------------------------------------------------------//

process extractNonUMIreads {
  tag "${prefix}"
  label 'samtools'
  label 'medCpu'
  label 'extraMem'

  publishDir "${params.outDir}/extractNonUMIreads", mode: 'copy'

  input :
  set val(prefix), file(sortedBam), file(nonUmisReadsIDs) from chAlignBamCheck_allreads.join(chUmiReadsIDs_exNonUMIreads)

  output:
  set val("${prefix}"), file("*_NonUmi_assignedNonUMIs.{bam,bam.bai}") into chNonUmiBam, chNonUmiBam_geneBody

  script:  
  """
  # Separate umi and non umi reads
  samtools view ${sortedBam[0]} > ${prefix}assignedAll.sam

  # save header and extract non umi reads 
  samtools view -H ${sortedBam[0]} > ${prefix}_assignedNonUMIs.sam

  nbLines=\$(wc -l < ${nonUmisReadsIDs})
  # get reads that match non umi read IDs
  if((\$nbLines!=0))
  then
    fgrep -f ${nonUmisReadsIDs} ${prefix}assignedAll.sam >> ${prefix}_assignedNonUMIs.sam || echo "no sequence found"
  fi
 
  # sam to bam
  samtools view -bh ${prefix}_assignedNonUMIs.sam > ${prefix}_NonUmi_assignedNonUMIs.bam

  # index
  samtools index ${prefix}_NonUmi_assignedNonUMIs.bam
  rm *.sam
  """
}

// All reads  -------------------------------------------------------------------------------------------------------//

// Remove PCR duplicates  -------------------------------------------------------------------------------//

process chRmPcrDup_samtools {
  tag "${prefix}" // "${prefix}_NonUmi"
  label 'samtools'
  label 'medCpu'
  label 'medMem'

  publishDir "${params.outDir}/rmPcrDup", mode: 'copy'

  input:
  set val(prefix), file(aligned) from chNonUmiBam

  output:
  set val(prefix), file("*_dedup.bam") into chNonUmi_dedup
  set val(prefix), file("*_dedup.log") into chNonUmi_dedup_log
  set val(prefix),file("*dedup_summary.log") into chNonUmi_dedup_mqc

  script :
  """
    samtools collate ${aligned[0]} -o namecollate.bam 
    ##Add ms and MC tags for markdup to use later:
    samtools fixmate -m namecollate.bam  fixmate.bam
    #Markdup needs position order:
    samtools sort fixmate.bam -o positionsort.bam 
    #Finally mark duplicates:
    samtools markdup positionsort.bam ${prefix}_NonUmi_dedup.bam -s -r &> ${prefix}_NonUmi_dedup.log
    
    rm fixmate.bam positionsort.bam namecollate.bam

    # get percent dup
    dup=\$(grep "DUPLICATE TOTAL:"  ${prefix}_NonUmi_dedup.log | cut -f3 -d" ") 
    tot=\$(grep "READ:" ${prefix}_NonUmi_dedup.log | cut -f2 -d" ")
    percent_dup=\$(echo "\$dup \$tot" | awk ' { printf "%.*f", 2, \$1/\$2 } ')
    # pour plot dans mqc : prefix, x, y
    # x=number of duplicates, y=percent of duplicates
    echo ${prefix} "," \$tot "," \$dup "," \$percent_dup > ${prefix}_NonUmi_dedup_summary.log
  """
}

process chRmPcrDup_umitools {
  tag "${prefix}" // "${prefix}_Umi"
  label 'umiTools'
  label 'medCpu'
  label 'medMem'

  publishDir "${params.outDir}/rmPcrDup", mode: 'copy'

  input:
  set val(prefix), file(umiBam) from chAssignBam_dedup

  output:
  set val(prefix), file("*_dedup.bam") into chUmi_dedup
  set val(prefix), file("*_dedup.log") into chUmi_dedup_log
  set val(prefix), file("*_dedup_summary.log") into chUmi_dedup_mqc


  script:
  """
  # Dedup per position and UMI (hamming=1)
  umi_tools dedup --stdin=${umiBam[0]}  --log=${prefix}_Umi_dedup.log > ${prefix}_Umi_dedup.bam


  # get percent dup
    tot=\$(grep "Reads: Input Reads:"  ${prefix}_Umi_dedup.log | cut -f7 -d" " | sed 's/,//') 
    dedup=\$(grep "Number of reads out:" ${prefix}_Umi_dedup.log | cut -f8 -d" " | sed 's/,//') 
    percent_dedup=\$(echo "\$dedup \$tot" | awk ' { printf "%.*f", 2, \$1/\$2*100 } ')
    # pour plot dans mqc : prefix, x, y
    # x=number of dedduplicates, y=percent of duplicates
    echo ${prefix}","\$tot","\$dedup","\$percent_dedup > ${prefix}_Umi_dedup_summary.log
  """
}


// Merge umi & non umi  -----------------------------------------------------------------//


process chMergeUmiNonUmiBam {
  tag "${prefix}"
  label 'samtools'
  label 'highCpu'
  label 'medMem'

  publishDir "${params.outDir}/featurecounts/allreads", mode: 'copy'

  input:
  set val(prefix), file(umibam), file(nonumibam) from chUmi_dedup.join(chNonUmi_dedup)

  output:
  set val(prefix), file("*_merged_dedup.bam") into chMergeDedupBam
 
  script:
  """
  samtools merge -@ ${task.cpu} ${prefix}_merged_dedup.bam ${umibam} ${nonumibam}
  """
}

// Assign  -------------------------------------------------------------------------------//

process readAssignmentAllReads {
  tag "${prefix}"
  label 'featureCounts'
  label 'highCpu'
  label 'medMem'

  publishDir "${params.outDir}/featurecounts/allreads", mode: 'copy'

  input :
  set val(prefix), file(bam) from chMergeDedupBam
  file(genome) from chGtfFCallreads.collect()

  output : 
  set val(prefix), file("*featureCounts.bam") into chAssignMergeDedupBam
  file "*.summary" into chAssignmentLogs_allreads
  set val(prefix), file("*_counts") into featureCountMatrix_allreads

  script:
  """	
  featureCounts  -p \
    -a ${genome} \
    -o ${prefix}_counts \
    -T ${task.cpus} \
    -R BAM \
    -g gene_name ${bam}

  featureCounts -v &> v_featurecounts.txt

  # -a annotation file
  # -R results format
  """
}

process sortAndIndexBamAllReads {
  tag "${prefix}"
  label 'samtools'
  label 'highCpu'
  label 'medMem'

  publishDir "${params.outDir}/bams/allreads", mode: 'copy'

  input:
  set val(prefix), file(assignDedupBam) from chAssignMergeDedupBam
	
  output:
  set val(prefix), file("*_Sorted.{bam,bam.bai}") into chSortedBAMBigWig, chSortedBAMSaturationCurve

  script :
  """
  samtools sort -@ ${task.cpus} ${assignDedupBam} -o ${prefix}_Sorted.bam
  samtools index ${prefix}_Sorted.bam
  """
}

//percent dup = (dup seq umi + dup seq non umi)/tot seq before dedup

// Matrix  -------------------------------------------------------------------------------//

// per cell
process countMatricesAllReads {
  tag "${prefix}"
  label 'featureCounts'
  label 'lowCpu'
  label 'highMem'

  input:
  set val(prefix), file(featureCountsBed) from featureCountMatrix_allreads

  output:
    set val(prefix), file("*_readCounts.tsv.gz") into chMatricesRead
    set val(prefix), file("*_nbGenes.txt") into chReadCountGenes //for mqc 

  script:
  """
    grep -v "^#"  ${featureCountsBed} | cut -f 1,7 | tail -n+2 >> ${prefix}"_selected"
    awk '{if(\$2!=0) print }' ${prefix}"_selected" >> ${prefix}"_readCounts.tsv"
    nbgene=\$(wc -l ${prefix}"_readCounts.tsv" | cut -f1 -d" ")
    echo ${prefix} \$nbgene > ${prefix}"_nbGenes.txt"
    echo -e 'gene count' > header
    cat ${prefix}"_readCounts.tsv" >> header
    mv header ${prefix}"_readCounts.tsv"
    gzip ${prefix}"_readCounts.tsv"
  """
}

// merge all cell
process mergeReadMatrices {
  tag "${prefix}"
  label 'R'
  label 'highCpu'
  label 'medMem'

  publishDir "${params.outDir}/Matrices", mode: 'copy'

  input:
  file (readmatrices) from chMatricesRead.collect()

  output:
  file ("10XlikeMatrix_read") into ch10X_allreads
  file ("10XlikeMatrix_read.zip") into ch10Xzip_allreads
  file ("read_gene_per_cell.txt") into chReadResume

  script:
  """
  merge_matrices.r 10XlikeMatrix_read/ read
  zip 10XlikeMatrix_read.zip 10XlikeMatrix_read/*
  """ 
}

// Bigwig  ----------------------------------------------------------------------//

process bigWig {
  tag "${prefix}"
  label 'deeptools'
  label 'extraCpu'
  label 'medMem'

  publishDir "${params.outDir}/AllReads/bigWig", mode: 'copy'

  input:
  set val(prefix), file(bam) from chSortedBAMBigWig

  output:
  set val(prefix), file("*_coverage.bw") into chBigWig 
  set val(prefix), file("*_coverage.log") into chBigWigLog
  file("v_deeptools.txt") into chBamCoverageVersion

  script:
  """
  ## Create bigWig files
  bamCoverage --normalizeUsing CPM -b ${bam[0]} -of bigwig -o ${prefix}_coverage.bw --numberOfProcessors=${task.cpus}  > ${prefix}_coverage.log
  bamCoverage --version &> v_deeptools.txt
  """
}

// Analysis all reads ----------------------------------------------------------------------//

process saturationCurves {
  tag "${prefix}"
  label 'preseq'
  label 'extraCpu'
  label 'extraMem'

  errorStrategy 'ignore'

  when:
  !params.skipSatCurves

  input:
  set val(prefix), file(sortBam) from chSortedBAMSaturationCurve

  output:
  set val(prefix), file ("*curve.txt") into preseq_results
  file("v_preseq.txt") into chPreseqVersion

  script:
  """
  preseq lc_extrap -v -B ${sortBam[0]} -o ${prefix}.extrap_curve.txt -e 200e+06 &> ${prefix}.extrap_curve.log

  if grep ERROR ${prefix}.extrap_curve.log
  then 
    touch ${prefix}.extrap_curve.txt
  fi
  
  # install bedtools
  # bedtools bamtobed [OPTIONS] -i <BAM>
  # prendre le bed en input pour gc_extrap

  preseq gc_extrap 
  # -e, -extrap = Max extrapolation. Here extrapolate until 200 000 000 reads
  # -D, -defects = estimates the complexity curve without checking for instabilities in the curve.
  # -s, -step The step size for samples. Default is 1 000 000 reads
  # -n, -bootstraps The number of bootstraps. Default is 100
  # -c, -cval Level for confidence intervals. Default is 0.95
  # -d, -dupllevelFraction of duplicates to predict. Default is 0.5
  # -x, -termsMax number of terms for extrapolation. Default is 100
  # -pe,  Input is a paired end read file
  preseq &> v_preseq.txt
  """
}

// Gene-based saturation
process geneSaturation {
  label 'R'
  label 'medCpu'
  label 'medMem'

  //when:
  //!params.skip_qc && !params.skip_saturation

  input:
  file (matDir) from ch10X_allreads

  output:
  file "*gcurve.txt" into genesat_results

  script:
  """
  gene_saturation.r ${matDir} counts.gcurve.txt
  """
}

//Gene body Coverage
process genebodyCoverage {
  tag "${prefix}"
  label 'rseqc'
  label 'medCpu'
  label 'medMem'
  errorStrategy 'ignore'
  publishDir "${params.outDir}/genebody_coverage" , mode: 'copy',
  saveAs: {filename ->
      if (filename.indexOf("geneBodyCoverage.curves.pdf") > 0)       "geneBodyCoverage/$filename"
      else if (filename.indexOf("geneBodyCoverage.r") > 0)           "geneBodyCoverage/rscripts/$filename"
      else if (filename.indexOf("geneBodyCoverage.txt") > 0)         "geneBodyCoverage/data/$filename"
      else if (filename.indexOf("log.txt") > -1) false
      else filename
  }

  when:
  !params.skipGeneCov

  input:
  file bed12 from chBedGeneCov.collect()
  set val(prefix), file(bm) from chUmiBam_geneBody.concat(chNonUmiBam_geneBody) 

  output:
  file "*.{txt,pdf,r}" into chGeneCov_res
  file ("v_rseqc") into chRseqcVersion

  script:
  """
  geneBody_coverage.py \\
      -i ${bm[0]} \\
      -o ${prefix}.rseqc \\
      -r $bed12
  mv log.txt ${prefix}.rseqc.log.txt
  geneBody_coverage.py --version &> v_rseqc
  """
}

// UMI reads only --------------------------------------------------------------------------


/*process DupPerGeneRatio{
  tag "${prefix}"
  label 'R'
  label 'lowCpu'
  label 'lowMem'

  input:

  output:
  file ("DupPerGene.csv") into chDupPerGene_mqc

  script:
  """
  
  """ 
}
*/


process umiPerGeneDist{
  tag "${prefix}"
  label 'R'
  label 'lowCpu'
  label 'lowMem'

  publishDir "${params.outDir}/umiPerGeneDist", mode: 'copy'

  input:
  set val(prefix), file(matrix) from chMatrices_dist

  output:
  set val(prefix), file ("*_HistUMIperGene_mqc.csv") into chUMIperGene

  script:
  """
  # Get matrix one by one
  umiPerGene_dist.r ${matrix} ${prefix}
  """ 
}

process mtRNAratio {
  tag "${prefix}"
  label 'R'
  label 'highCpu'
  label 'medMem'

  input:
  file (matDir) from ch10X_mt
  file (umiSummary) from chUmiResume_mt

  output:
  file ("RatioPerCell.csv") into chUmiGeneRatio
  file ("MtGenePerCell.csv") into chMT

  script:
  """
  mt_ratio_rna.r ${matDir} ${params.genome} ${umiSummary}
  """ 
}


process countUMIGenePerCell{
  tag "${prefix}"
  label 'R'
  label 'lowCpu'
  label 'medMem'

  input:
  file(matrices) from chMatrices_UMIGenePerCell.collect()

  output:
  file ("nbGenePerCell.csv") into chGenePerCell
  file ("nbUMIPerCell.csv") into chUmiPerCell

  script:
  """
  umiGenePerCell.r
  """ 
} 

// MultiQC ----------------------------------------------------------------------//

process getSoftwareVersions{
  label 'python'
  label 'lowCpu'
  label 'lowMem'
  publishDir path: "${params.outDir}/softwareVersions", mode: "copy"

  when:
  !params.skipSoftVersions

  input:
  file("v_umi_tools.txt") from chUmiToolsVersion.first().ifEmpty([])
  file("v_cutadapt.txt") from chCutadaptVersion.first().ifEmpty([])
  file("v_star.txt") from chStarVersion.first().ifEmpty([])
  file("v_featurecounts.txt") from chFCversion.first().ifEmpty([])
  file("v_samtools.txt") from chSamtoolsVersion.first().ifEmpty([])
  file("v_deeptools.txt") from chBamCoverageVersion.first().ifEmpty([])
  file("v_R.txt") from chRversion.ifEmpty([])
  file("v_rseqc") from chRseqcVersion.first().ifEmpty([])
  file("v_preseq.txt") from chPreseqVersion.first().ifEmpty([])

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
  label 'lowCpu'
  label 'lowMem'
  label 'multiqc'

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
  section_href: 'https://gitlab.curie.fr/sc-platform/smartseq3'
  plot_type: 'html'
  data: |
        <dl class=\"dl-horizontal\">
  ${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
  """.stripIndent()
}

if (skippedPoorAlignment.size() > 0){
  Channel.fromList(skippedPoorAlignment)
         .flatMap{ it -> it + ": Poor alignment rate. Sample removed from the analysis !!!! <br>"}
         .collectFile(name: 'warnings.txt', newLine: true)
         .set{chWarn}
}else{
  chWarn = Channel.empty()
}

process multiqc {
  label 'multiqc'
  label 'medCpu'
  label 'lowMem'
  publishDir "${params.outDir}/MultiQC", mode: 'copy'

  when:
  !params.skipMultiQC

  input:
  file splan from chSplan.collect()
  file multiqcConfig from chMultiqcConfig
  file ('software_versions/*') from softwareVersionsYaml.collect().ifEmpty([])
  file ('workflowSummary/*') from workflowSummaryYaml.collect()
  file ('workflowSummary/*') from chWarn.collect().ifEmpty([]) 
  //Modules
  file ('star/*') from chAlignmentLogs.collect().ifEmpty([])
  file ('FC/*') from chAssignmentLogs_allreads.collect().ifEmpty([])
  file ('coverage/*') from chGeneCov_res.collect().ifEmpty([])
  file ('preseq/*') from preseq_results.collect().ifEmpty([])
  //LOGS
  //file ('umiExtract/*') from chUmiExtractedLog.collect()
  file('pUMIs/*') from chCountSummaryExtUMI.collect()
  file('totFrag/*') from chTotFrag.collect()
  file ('bigwig/*') from chBigWigLog.collect()
  file (umiResume) from chUmiResume // general stats 
  file (readResume) from chReadResume
  //PLOTS
  file ("umiPerGene/*") from chUMIperGene.collect() // linegraph == histogram
  file ("nbUMI/*") from chUmiPerCell.collect()  // bargraph
  file ("nbGene/*") from chGenePerCell.collect() // bargraph 
  file ("ratio/*") from chUmiGeneRatio.collect() // dotplot UmiGenePerCell_mqc.csv
  file ("mt/*") from chMT.collect() // dotplot MtGenePerCell_mqc.csv
  //file ("dupPerGene/*") from chDupPerGene_mqc.collect() // dotplot DupPerGene.csv
  file ('genesat/*') from genesat_results.collect().ifEmpty([]) // linegraph
  file('rmPCRumi/*') from chUmi_dedup_mqc.collect()


  output: 
  file splan
  file "*report.html" into multiqc_report
  file "*_data"

  script:
  rtitle = customRunName ? "--title \"$customRunName\"" : ''
  rfilename = customRunName ? "--filename " + customRunName + "_report" : "--filename report"
  modules_list = "-m custom_content -m samtools -m star -m featureCounts -m deeptools -m preseq -m rseqc"
  //warn=skippedPoorAlignment.size() > 0 ? "--warn workflowSummary/warnings.txt" : ""
  """
  stat2mqc.sh ${splan}
  #mean_calculation.r
  mqc_header.py --splan ${splan} --name "${params.protocol} scRNA-seq" --version ${workflow.manifest.version} > multiqc-config-header.yaml
  multiqc . -f $rtitle $rfilename -c multiqc-config-header.yaml -c $multiqcConfig $modules_list --cl_config '{max_table_rows: 2000}'
  """
}


//Sub-routines
process outputDocumentation {
  label 'python'
  label 'lowCpu'
  label 'lowMem'

  publishDir "${params.outDir}/summary", mode: 'copy'

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
  reportFields['skippedPoorAlignment'] = skippedPoorAlignment
  
  // Render the TXT template
  def engine = new groovy.text.GStringTemplateEngine()
  def tf = new File("$projectDir/assets/onCompleteTemplate.txt")
  def txtTemplate = engine.createTemplate(tf).make(reportFields)
  def reportTxt = txtTemplate.toString()

  // Render the HTML template
  def hf = new File("$projectDir/assets/onCompleteTemplate.html")
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
  File woc = new File("${params.outDir}/workflowOnComplete.txt")
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

  //final logs
  if(skippedPoorAlignment.size() > 0){
    log.info "WARNING - ${skippedPoorAlignment.size()} samples skipped due to poor alignment scores!"
  }

  if(workflow.success){
      log.info "Pipeline Complete"
  }else{
    log.info "FAILED: $workflow.runName"
  }
}
