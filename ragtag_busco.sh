#!/usr/bin/env bash

################################################################################
# Help                                                                         #
################################################################################
Help()
{
   # Display Help
   usage="
   ragtag_busco automates the 3 steps of ragtag (correct, scaffold, patch)
   to scaffold a draft genome based on a reference then runs a BUSCO on the
   result draft.
   Syntax: ragtag_busco -g draft_genome -r reference_genome -o output_dir
   options:
   -h     Print this Help.
   -r     Reference file to scaffold on.
   -g     Genome file to scaffold.
   -o     Output directory.
   -l     Lineage for BUSCO (default : eukaryota)."
}

################################################################################
################################################################################
# Main program                                                                 #
################################################################################
################################################################################

Help
echo "Ragtag_busco"

lineage="eukaryota"

while getopts :h::r:g:o:l: flag
do
    case "${flag}" in
        h) echo "$usage"; exit;;
        r) reference=${OPTARG};;
        g) genome=${OPTARG};;
        o) output=${OPTARG};;
        l) lineage=${OPTARG};;
        :) printf "missing argument for -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
       \?) printf "illegal option: -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
    esac
done

# mandatory arguments
if [ ! "$reference" ] || [ ! "$genome" ] || [ ! "$output" ]; then
  echo "arguments -r, -g and -o must be provided"
  echo "$usage" >&2; exit 1
fi

echo "######"
echo "### Ragtag step 1/3 : Correct"
echo "######"

ragtag.py correct ${reference} ${genome} -o ${output}/correct

echo "######"
echo "### Ragtag step 2/3 : scaffold"
echo "######"
ragtag.py scaffold ${reference} ${output}/correct/ragtag.correct.fasta -r -o ${output}/scaffold

echo "######"
echo "### Ragtag step 3/3 : patch"
echo "######"
ragtag.py patch ${reference} ${output}/scaffold/ragtag.scaffold.fasta -o ${output}/patch

echo "######"
echo "### BUSCO"
echo "######"
busco -i ${output}/patch/ragtag.patch.fasta -o busco --mode genome --lineage_dataset ${lineage} --tar -f
