#!/usr/bin/env bash

################################################################################
################################################################################
# transcriptomic assembler pipeline                                            #
################################################################################
################################################################################


################################################################################
# Functions                                                                    #
################################################################################
Help()
{
   # Display Help
   echo "
   #######
   ### Transcriptomic_pipeline
   #######

   Not done yet so here's some Dylan's lyrics waiting for a final help :

   There must be some way out of here
   Said the joker to the thief
   There's too much confusion
   I can't get no relief

   Syntax:
   transcriptomic_pipeline -M [C|R|T|A] -o output_dir -F forward_reads -R reverse_reads [-h|E|v|V|t|m|D|p|l|M|u]
   options:
   -h         Print this Help.
   ### Mandatory
   -M [C|R|T|A] working mode :
              C   check : check md5
              R   reads : check reads quality
              T   trimm : trim & filter then assess reads quality
              A   assemble : assemble and evaluate transcriptome
   -o dir     Output directory.
   ### Inputs
   -F reads1,reads2,reads3  Forward reads file (separated by comma)
   -R reads1,reads2,reads3  Reverse reads file (separated by comma)
   ### Exec options
   -t X       Number of threads
   -m X       Max memory (in G)
   -v         Verbose mode
   -V         Print version and exit
   ### Options
   -f         Force indexes/database generation (erase previous databases)
   -d         Download missing db
   -q         Run multiQC
   -u file    Fasta file of blascklisted RNA to remove
   -s file    Fasta file of a protein DB to BLAST against
   -l lineage Lineage for BUSCO (default : eukaryota).

   examples :
        heresanexample
        heresasecondexample
        heresathirdexample
   "
}

Version()
{
  echo "transcriptomic_pipeline v0.1"
}

Fastqc()
{
  # Exocet specific : cannot use conda's java #
  source deactivate
  #############################################

  # Reads quelity evaluation tool
  echo "starting fastqc"

  # Verbose
  if [[ ${verbose} ]]; then
    echo mkdir -p ${output}/fastqc_reports
    echo fastqc -t ${num_threads_fastqc} ${splittedFiles} -o ${output}/fastqc_reports/
  fi

  mkdir -p "${output}/fastqc_reports"
  # Fastqc allow a max of 6 threads
  if [[ ${num_threads} -gt 6 ]]; then
    num_threads_fastqc=6
  else
    num_threads_fastqc=${num_threads}
  fi
  fastqc ${splittedFiles} -o ${output}/fastqc_reports/

  # Exocet specific : cannot use conda's java #
  conda activate
  #############################################

  if [ "${run_multiqc}" ]; then
    echo "running multiqc"
    multiqc ${output}/fastqc_reports/ -o ${output}/fastqc_reports/
  fi
  echo "fastqc done"
}

Cutadapt()
{
  # Reads trimming tool : removing adapters
  echo "starting cutadapt"
  mkdir -p "${output}/cutadapt"
  cutadapt \
      -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
      -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
      -o ${output}/cutadapt/adapters.${name_forward} \
      -p ${output}/cutadapt/adapters.${name_reverse} \
      "${process_forward}" \
      "${process_reverse}" ;
  echo "cutadapt done"
}

Trimmomatic()
{
  # Reads trimming tool : quality trimming
  echo "Starting trimmomatic"
  mkdir -p "${output}/trimmomatic"
  trimmomatic PE \
      ${output}/cutadapt/adapters.${name_forward} \
      ${output}/cutadapt/adapters.${name_reverse} \
      ${output}/trimmomatic/paired_trimmed.${name_forward} \
      ${output}/trimmomatic/unpaired_trimmed.${name_forward} \
      ${output}/trimmomatic/paired_trimmed.${name_reverse} \
      ${output}/trimmomatic/unpaired_trimmed.${name_reverse} \
      LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 ;
  echo "Trimmomatic done"
}

Bowtie2_build()
{
  echo "Starting building Bowtie2 index"
  bowtie2-build ${file_to_index} ${file_to_index}
  echo "Bowtie2 index built"
}

Bowtie2_filter_align()
{
  echo "Starting bowtie2 alignement"
  mkdir -p "${output}/bowtie2"
  bowtie2 --quiet --very-sensitive-local --phred33 --threads ${num_threads} \
    -x "${unwantedRNA_file}" \
    -1 ${output}/trimmomatic/paired_trimmed.${name_forward} \
    -2 ${output}/trimmomatic/paired_trimmed.${name_reverse} \
    --met-file ${output}/bowtie2/${name_forward}_bowtie2Metrics.txt \
    --al-conc-gz ${output}/bowtie2/blacklist_paired_aligned.${name_forward} \
    --un-conc-gz ${output}/bowtie2/blacklist_paired_unaligned.${name_forward} \
    --al-gz ${output}/bowtie2/blacklist_unpaired_aligned.${name_forward} \
    --un-gz ${output}/bowtie2/blacklist_unpaired_unaligned.${name_forward} \

  mv ${output}/bowtie2/blacklist_paired_unaligned.${name_forward::-3}.1.gz\
    ${output}/bowtie2/filtered.${name_forward}
  mv ${output}/bowtie2/blacklist_paired_unaligned.${name_forward::-3}.2.gz\
    ${output}/bowtie2/filtered.${name_reverse}
  echo "Bowtie2 alignement done"
}

Trinity()
{
  echo "Starting Trinity assembler"
  mkdir -p "${output}/trinity"
  Trinity --seqType fq --max_memory ${max_memory} --CPU ${num_threads} \
          --output "${output}/trinity" \
          --left ${forward_file} \
          --right ${reverse_file}
  echo "Trinity assembly done"
}

Bowtie2_reads_representation()
{
  echo "Starting evaluating reads representation"
  mkdir -p ${output}/assembly_evaluation/reads_representation
  bowtie2 --threads ${num_threads} --phred33 --no-unal -k 20\
          -x ${output}/trinity/Trinity.fasta\
          -1 "${process_forward}"\
          -2 "${process_reverse}"\
          2>${output}/assembly_evaluation/reads_representation/${name_forward}_align_stats.txt|\
          samtools view -@10 -Sb -o ${output}/assembly_evaluation/reads_representation/${name_forward}_bowtie2.bam;

  tail -n 1 ${output}/assembly_evaluation/reads_representation/${name_forward}_align_stats.txt

  samtools sort ${output}/assembly_evaluation/reads_representation/${name_forward}_bowtie2.bam \
  -o ${output}/assembly_evaluation/reads_representation/${name_forward}_bowtie2.coordSorted.bam;
  samtools index ${output}/assembly_evaluation/reads_representation/${name_forward}_bowtie2.coordSorted.bam;

  echo "Evaluating reads representation done. You may want to open
  ${name_forward}_bowtie2.coordSorted.bam and Trinity.fasta with IGV to visualize
  alignement."
}

Blast_build()
{
  echo "Starting building Blast index"
  makeblastdb -in ${proteinDB_file} -dbtype prot
  echo "Blast index built"
}

Full_length_transcripts()
{
  echo "Starting Blast"
  mkdir -p ${output}/assembly_evaluation/blast
  blastx -evalue 1e-20 -num_threads ${num_threads} -max_target_seqs 1 -outfmt 6\
          -query ${output}/trinity/Trinity.fasta
          -db ${proteinDB_file}\
          -out ${output}/assembly_evaluation/blast/blastx.outfmt6
  echo "Blast done"
  echo "Starting trinity's analyze_blastPlus_topHit_coverage's script"
  $TRINITY_HOME=`which Trinity`
  $TRINITY_HOME/util/analyze_blastPlus_topHit_coverage.pl\
                  ${output}/assembly_evaluation/blast/blastx.outfmt6\
                  ${output}/trinity/Trinity.fasta\
                  ${proteinDB_file}
  echo "analyze_blastPlus_topHit_coverage done"
}

Busco()
{
  echo "Starting Busco"
  mkdir -p ${output}/assembly_evaluation/busco
  busco\
        -i ${output}/trinity/Trinity.fasta\
        -o ${output}/assembly_evaluation/busco --mode transcriptome\
        --lineage_dataset ${lineage} --tar --download_path ${output}/reference/
  echo "Busco done"
}

################################################################################
# Main program                                                                 #
################################################################################

### Argument parsing

while getopts hM:o:F:R:t:m:vVfdql:u:s: options
do
    case "${options}" in
        h) Help; exit;;
        V) Version; exit;;
        M) mode=${OPTARG};;
        o) output=${OPTARG};;
        F) forward_file=${OPTARG};;
        R) reverse_file=${OPTARG};;
        t) num_threads=${OPTARG};;
        m) max_memory=${OPTARG};;
        v) verbose=True;;
        f) force_indexes=True;;
        d) download_db=True;;
        q) run_multiqc=True;;
        l) lineage=${OPTARG};;
        u) unwantedRNA_file=${OPTARG}
           if [ ! -f "${unwantedRNA_file}" ]; then
             echo "${unwantedRNA_file} does not exists"; exit1;
           fi;;
        s) proteinDB_file=${OPTARG}
          if [ ! -f "${proteinDB_file}" ]; then
            echo "${proteinDB_file} does not exists"; exit1;
          fi;;

    esac
done

### Arguments checking
# mandatory arguments
if [ ! "${mode}" ] || [ ! "${output}" ] || [ ! "${forward_file}" ] || [ ! "${reverse_file}" ]; then
  echo "arguments -M, -o, -F and -R must be provided"
  Help >&2; exit 1
fi

# Default for threads and memory
if [[ ! ${num_threads} ]]; then
  num_threads=1
fi

if [[ ! ${max_memory} ]]; then
  max_memory=16
fi

if [[ ! ${lineage} ]]; then
  lineage="eukaryota"
fi

# Merge all input files in 1 variable
splittedFiles="$(echo ${forward_file}|tr ',' ' ') $(echo ${reverse_file}|tr ',' ' ')"
# Check if input files exists
for input_file in ${splittedFiles}; do
  if [ ! -f "${input_file}" ]; then
    echo "${input_file} does not exists. Please check your inputs"
    exit 1
  fi
done

# output dir in right format
output=$(echo "${output}" | sed 's:/$::')
output=$(echo "${output}" | sed 's: :_:')


case "${mode}" in

  ### MODE CHECK MD5
  C)  echo "a ajouter ..."; exit 0;;

  ### MODE CHECK READS QUALITY
  R)  echo "Starting in R mode : check reads quality"

      # Run Fastqc
      Fastqc

      # Done
      echo "R mode over. Check your fastqc reports files. You can then move on
to T mode."
      exit 0;;

  ### MODE TRIM & FILTER
  T)  echo "Starting in T mode : trim & filter then assess reads quality"

      # Dir creation for the reference files
      mkdir -p "${output}/reference"

      # Get input files in an array to be able to parse them
      IFS=',' read -ra forwardFiles_array <<< ${forward_file}
      IFS=',' read -ra reverseFiles_array <<< ${reverse_file}

      # Verification of the number of input files
      if [[ ${#forwardFiles_array[*]} -ne ${#reverseFiles_array[*]} ]]; then
        echo "You have a different number of forward and reverse files. Please
        check your files"; Help >&2; exit 1;
      fi

      # Verification of the unwantedRNA_file
      if [[ ! "${unwantedRNA_file}" ]]; then

        # Won't proceed without an unwantedRNA_file
        if [[ ! "${download_db}" ]] && [[ ! -f ${output}/reference/unwanted_RNA.fa.gz ]]; then
          echo "You need to provide a file of RNA to filter out with -u file or
          allow the programm to download one with -d"; Help >&2; exit 1;

        # Download of the db if asked (2 files : LSU and SSU + 2 md5 files to
        # check download integrity)

        # SSU
        elif [[ ! -f ${output}/reference/unwanted_RNA.fa.gz ]]; then
          wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/SILVA_138.1_SSURef_tax_silva_trunc.fasta.gz \
          -P ${output}/reference/
          wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/SILVA_138.1_SSURef_tax_silva_trunc.fasta.gz.md5 \
          -P ${output}/reference/
          # Verification of the integrity of the downloaded file
          if [[ $( md5sum ${output}/reference/SILVA_138.1_SSURef_tax_silva_trunc.fasta.gz|cut -f 1 -d ' ') != $(cat ${output}/reference/SILVA_138.1_SSURef_tax_silva_trunc.fasta.gz.md5|cut -f 1 -d ' ') ]]; then
            echo "There was a problem with the download. Please try again or use
            your own file with -u file."; exit 1;
          fi

          # LSU
          wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/SILVA_138.1_LSURef_tax_silva_trunc.fasta.gz \
          -P ${output}/reference/
          wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/SILVA_138.1_LSURef_tax_silva_trunc.fasta.gz.md5 \
          -P ${output}/reference/
          # Verification of the integrity of the downloaded file
          if [[ $( md5sum ${output}/reference/SILVA_138.1_LSURef_tax_silva_trunc.fasta.gz|cut -f 1 -d ' ') != $(cat ${output}/reference/SILVA_138.1_LSURef_tax_silva_trunc.fasta.gz.md5|cut -f 1 -d ' ') ]]; then
            echo "There was a problem with the download. Please try again or use
            your own file with -u file."; exit 1;
          fi

          # Merging both SILVA files into 1 unwantedRNA_file
          # SILVA database contains U instead of T. Will modify the header
          # aswell but it doesn't matter in our case
          echo "Creation of unwanted_RNA.fa.gz"
          touch ${output}/reference/unwanted_RNA.fa
          zcat ${output}/reference/SILVA_138.1_LSURef_tax_silva_trunc.fasta.gz|sed 's/U/T/g' > ${output}/reference/unwanted_RNA.fa
          zcat ${output}/reference/SILVA_138.1_SSURef_tax_silva_trunc.fasta.gz|sed 's/U/T/g' >> ${output}/reference/unwanted_RNA.fa
          gzip -v ${output}/reference/unwanted_RNA.fa
          unwantedRNA_file="${output}/reference/unwanted_RNA.fa.gz"

        else
          echo "using found ${output}/reference/unwanted_RNA.fa.gz file"
          unwantedRNA_file="${output}/reference/unwanted_RNA.fa.gz"
        fi

      fi

      # Building Bowtie2 unwantedRNA index if -f or none found
      if [[ ${force_indexes} ]] || [[ ! -f ${unwantedRNA_file}.1.bt2 ]]; then
        file_to_index=${unwantedRNA_file}
        Bowtie2_build
      else
        echo "Using found bowtie2 index"
      fi

      splittedFiles=""
      for (( n=0; n < ${#forwardFiles_array[*]}; n++ )); do
        echo "Processing pair or reads number ${n}"
        process_forward=${forwardFiles_array[n]}; name_forward="$(basename -- ${process_forward})"
        process_reverse=${reverseFiles_array[n]}; name_reverse="$(basename -- ${process_reverse})"
        Cutadapt
        Trimmomatic
        Bowtie2_filter_align
        splittedFiles="${splittedFiles} ${output}/bowtie2/filtered.${name_forward} ${output}/bowtie2/filtered.${name_reverse}"
      done

      # Quality eval of the new reads with fastqc
      Fastqc

      echo "T mode is over. Check your new fastqc reports files. If you're happy
with your reads quality, move on to A mode."
      exit 0;;

  ### MODE CHECK MD5
  A)  echo "Starting in A mode : assemble and evaluate transcriptome"

      # Running Trinity assembler
      Trinity

      # Quality evaluation
      mkdir -p ${output}/assembly_evaluation

      # Get input files in an array to be able to parse them
      IFS=',' read -ra forwardFiles_array <<< ${forward_file}
      IFS=',' read -ra reverseFiles_array <<< ${reverse_file}

      # Verification of the number of input files
      if [[ ${#forwardFiles_array[*]} -ne ${#reverseFiles_array[*]} ]]; then
        echo "You have a different number of forward and reverse files. Please
        check your files"; Help >&2; exit 1;
      fi

      # Reads representation
      ${file_to_index}="${output}/trinity/Trinity.fasta"
      Bowtie2_build
      for (( n=0; n < ${#forwardFiles_array[*]}; n++ )); do
        echo "Processing pair or reads number ${n}"
        process_forward=${forwardFiles_array[n]}; name_forward="$(basename -- ${process_forward})"
        process_reverse=${reverseFiles_array[n]}; name_reverse="$(basename -- ${process_reverse})"
        Bowtie2_reads_representation
      done

      # Counting Full Length Trinity Transcripts
      # Won't proceed without a uniprot/swissprot db
      if [[ ! "${proteinDB_file}" ]]; then
        if [[ ! "${download_db}" ]] && [[ ! -f ${output}/reference/protein_DB.fa.gz ]]; then
          echo "No protein database found. Please provide a protein db file to
          blast against with -s file or allow download with -d"; Help >&2; exit 1;

      elif [[ ! -f ${output}/reference/protein_DB.fa.gz ]]; then
          wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz \
          -o ${output}/reference/protein_DB.fa.gz
      fi
      proteinDB_file=${output}/reference/protein_DB.fa.gz
    fi

    if [[ ${force_indexes} ]] || [[ ! -f ${proteinDB_file::-3}.pdb ]]; then
      Blast_build
    fi

    Full_length_transcripts

    # Evaluating Busco score
    Busco

    echo "A mode is over. Your assembly is in ${output}/trinity/Trinity.fasta"
    exit 0;;

  *) echo "invalid mode. Please entrer a value in C,A,R,T"; Help >&2; exit 1;;
esac
#busco -i ${output}/patch/ragtag.patch.fasta -o busco --mode genome --lineage_dataset ${lineage} --tar -f
