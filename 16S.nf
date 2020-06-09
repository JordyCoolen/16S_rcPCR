params.threads = 4
params.outDir = "./results"
params.reads = "$baseDir/data/reads/test/ERR1664619_R{1,2}.fastq.gz"
params.database = "$baseDir/kraken2_db/minikraken_8GB_20200312/"
//params.KMA_database = "$baseDir/kma_db/NCBI_16S"
params.KMA_database = "$baseDir/kma_db/Demo_db2_k21"
params.blastdb_path = "$baseDir/blast_db/"
params.trim = true
params.krono = false
params.centrifuge = false
params.abricate = true

// Parsing the input parameters
sampleName       = "$params.sampleName"
outDir           = "$params.outDir" + "/" + "$params.sampleName"
threads          = "$params.threads"
trim             = "$params.trim"
krono            = "$params.krono"
centrifuge       = "$params.centrifuge"
abricate         = "$params.abricate"
awk              = "$baseDir/awk.bash"

// databases
database         = "$params.database"
KMA_database     = "$params.KMA_database"
blastdb_path     = "$params.blastdb_path"
centrifuge_database  = "$baseDir/centrifuge_db/p_compressed"

// special channel for the fastQ reads
Channel
      .fromFilePairs( params.reads )
      .ifEmpty { "cannot find read pairs in path"}
      .set  { reads_ch1 }

log.info """
NEXTFLOW 16S V0.50
================================
sample     : $params.sampleName
reads      : $params.reads
outDir     : $params.outDir
outFolder  : $outDir
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
threads    : $params.threads
database   : $database
KMA_db     : $KMA_database
blastdb    : $blastdb_path
trim       : $trim
fast       : $krono
centrifuge : $centrifuge
abricate   : $abricate
================================
"""

// Clean reads (adapter and read length filter)
process '1A_clean_reads' {
 tag '1A'
 conda 'bioconda::fastp=0.20.0 bioconda::fastqc=0.11.9 bioconda::prinseq=0.20.4'
 publishDir outDir + '/fastp', mode: 'copy'
 input:
   set pairID, file(reads) from reads_ch1
 output:
   //set file("${reads[0].baseName}_fastp.fastq.gz"), file("${reads[1].baseName}_fastp.fastq.gz") into fastp_1B, fastp_1D, fastp_1E, fastp_1F, fastp_1G, fastp_2
   set file("${sampleName}_R1.fastq.gz"), file("${sampleName}_R2.fastq.gz") into fastp_1B, fastp_1D, fastp_1E, fastp_1F, fastp_1G, fastp_2
   file "*"
   //file "${reads[0].baseName}_fastp.json" into fastp_json
   //file ".command.*"
 script:
   if( trim == true )
   """
   mkdir -p ${baseDir}/${outDir}/fastp/
   fastp -i ${reads[0]} -I ${reads[1]} -o ${reads[0].baseName}_fastp.fastq.gz -O ${reads[1].baseName}_fastp.fastq.gz \
   --length_required 50 --json ${reads[0].baseName}_fastp.json --html ${reads[0].baseName}_fastp.html --thread ${threads}
   fastqc -o ${baseDir}/${outDir}/fastp/ ${reads[0].baseName}_fastp.fastq.gz ${reads[1].baseName}_fastp.fastq.gz

   gunzip ${reads[0].baseName}_fastp.fastq.gz
   gunzip ${reads[1].baseName}_fastp.fastq.gz
   prinseq-lite.pl -derep 1 -derep_min 11 -fastq ${reads[0].baseName}_fastp.fastq -fastq2 ${reads[1].baseName}_fastp.fastq -out_good ${sampleName} -out_bad ${sampleName}_bad
   gzip ${reads[0].baseName}_fastp.fastq
   gzip ${reads[1].baseName}_fastp.fastq
   gzip ${sampleName}_bad_1.fastq
   gzip ${sampleName}_bad_2.fastq
   gzip ${sampleName}_1.fastq
   gzip ${sampleName}_2.fastq
   mv ${sampleName}_1.fastq.gz ${sampleName}_R1.fastq.gz
   mv ${sampleName}_2.fastq.gz ${sampleName}_R2.fastq.gz

   fastqc -o ${baseDir}/${outDir}/fastp/ ${sampleName}_R1.fastq.gz ${sampleName}_R2.fastq.gz

   """
   else if( trim == false )
   """
   mv ${reads[0]} ${sampleName}_R1.fastq.gz
   mv ${reads[1]} ${sampleName}_R2.fastq.gz
   mkdir -p ${baseDir}/${outDir}/fastp/
   fastqc -o ${baseDir}/${outDir}/fastp/ ${sampleName}_R1.fastq.gz ${sampleName}_R2.fastq.gz
"""
}

// Process 1B: Identify mycobacterium organisms (mainly for metagenomic data)
process '1B_ID_strain' {
  tag '1B'
  conda 'bioconda::kraken2=2.0.7_beta'
  publishDir outDir + '/kraken2', mode: 'copy'
  input:
  file reads from fastp_1B
  output:
    file "${sampleName}.kraken2.output.txt" into report_output
    file "${sampleName}.kraken2.report.txt" into report_file_1C, report_file_1D
    file ".command.*"
  script:
    """
    kraken2 --memory-mapping --confidence 0.05 --threads ${threads} --db ${database} --paired --report ${sampleName}.kraken2.report.txt \
    --output ${sampleName}.kraken2.output.txt ${reads[0]} ${reads[1]}
    """
}

// Process 1C: Krona plot
process '1C_Krona_plot' {
  tag '1B'
  conda 'bioconda::krona=2.7.1 anaconda::curl anaconda::make'
  publishDir outDir + '/krona', mode: 'copy'
  input:
  file output from report_output
  file report from report_file_1C
  output:
    file "*.html"
    file ".command.*"
  script:
    if ( krono == true )
    """
    ktUpdateTaxonomy.sh
    bash $awk ${report} > ${report}.edit
    ktImportTaxonomy ${report}.edit -o final_report.html
    """
    else if ( krono == false )
    """
    echo 'none' > final_report.html
    """
}

// Process 1D: Bracken
process '1D_bracken' {
  tag '1D'
  conda 'bioconda::bracken=2.5'
  publishDir outDir + '/kraken2', mode: 'copy'
  input:
  file report from report_file_1D
  output:
    file "${sampleName}.bracken.report.txt"
    file ".command.*"
  script:
    """
    bracken -d ${database} -i ${report} -o ${sampleName}.bracken.report.txt -l S
    """
}

// create KMA tool to detect 16S
// Process 1D: KMA
process '1D_KMA' {
  tag '1D'
  time "10m"
  conda 'bioconda::kma=1.2.22'
  publishDir outDir + '/kma', mode: 'copy'
  input:
  file reads from fastp_1D
  output:
    file "kma*"
    file "kma.fsa" into consensus_1F
    file "*.sam" into kma_1D
    file ".command.*"
  script:
    """
    #kma -t_db ${KMA_database} -mem_mode -ef -ex_mode -pm f -fpm f -sam 4 -ipe ${reads[0]} ${reads[1]} -o kma > kma.sam 2>/dev/null || exit 0
    # used for the 20200414_NimaGen_mix5_dilution_data, is more strict in selecting the correct reference which I think is best method
    kma -t_db ${KMA_database} -mem_mode -ef -ex_mode -ConClave 2 -apm f -sam 4 -ipe ${reads[0]} ${reads[1]} -o kma > kma.sam 2>/dev/null || exit 0
    """
}

// create KMA tool to match primers
// Process 1F: KMA
process '1F_process_KMA' {
  tag '1F'
  conda 'bioconda::samtools=1.9'
  errorStrategy 'ignore'
  publishDir outDir + '/kma', mode: 'copy'
  input:
  file sam from kma_1D
  output:
    file "${sampleName}.sorted.bam"
    file "${sampleName}.sorted.bam.bai"
    file ".command.*"
  script:
    """

    # merge sam files

    # sam --> bam
    samtools view -b ${sam} > ${sampleName}.bam
    # sort bam
    samtools sort ${sampleName}.bam > ${sampleName}.sorted.bam
    # index bam
    samtools index ${sampleName}.sorted.bam

    # additionally filter reference on only hits
    # this would make it possible to quickly evaluate the results

    """
}

// create abricate to detect 16S
// Process 1F: abricate
process '1F_abricate' {
  tag '1F'
  conda 'bioconda::abricate'
  publishDir outDir + '/abricate', mode: 'copy'
  input:
  file reads from fastp_1F
  file consensus from consensus_1F
  output:
    file "blast.txt"
    file ".command.*"
  script:
    if(abricate==true)
    """
    abricate --datadir ${blastdb_path} --db NCBI_16S ${consensus} --mincov 50 --minid 60  > blast.txt
    """
    else if(abricate==false)
    """
    echo 'none' > blast.txt
    """
}

// Process 1F: centrifuge
process '1G_centrifuge' {
  tag '1G'
  conda 'bioconda::centrifuge=1.0.4_beta'
  publishDir outDir + '/centrifuge', mode: 'copy'
  input:
  file reads from fastp_1G
  output:
    file "${sampleName}_report.txt" into report_file
    file ".command.*"
  script:
    if ( centrifuge == true )
    """
    centrifuge -p 4 -x ${centrifuge_database} --quiet -q -1 ${reads[0]} -2 ${reads[1]} \
    --report-file ${sampleName}_report.txt > ${sampleName}cent_stdout.out
    """
    else if ( centrifuge == false )
    """
    echo 'none' > ${sampleName}_report.txt
    """
}

// Process 2: multiQC
process '2_multiQC' {
  tag '2'
  publishDir outDir + '/QC', mode: 'copy'
  input:
  file reads from fastp_2
  output:
    file "*.html"
    file ".command.*"
  script:
    """
    multiqc ${baseDir}/${outDir}/fastp/
    """
}