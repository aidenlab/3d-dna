#!/bin/bash
## This is a script to set up initial conditions for replicating the AaegL4 assembly
## Input: None (link to NCBI assembly hardcoded)
## Output: AaegL2.fasta file
## Written by: OD
## USAGE: ./get-AaegL2.sh

[ -f temp.cprops ] || [ -f temp.asm ] && echo >&2 "Please remove or rename temp.cprops and/or temp.asm files from the current folder. Exiting!" && exit 1 

link="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/004/015/GCA_000004015.2_AaegL2/GCA_000004015.2_AaegL2_genomic.fna.gz"

pipeline=`cd "$( dirname $0)" && cd .. && pwd`

curl $link -o `basename $link`

curl `dirname $link`"/"`basename $link _genomic.fna.gz`"_assembly_report.txt" -o `basename $link _genomic.fna.gz`"_assembly_report.txt"

gunzip "`basename $link`"

awk -F "\t" '$1!~/#/&&$5!~/na/{count++; print $5, count, $9}' `basename $link _genomic.fna.gz`"_assembly_report.txt" > temp.cprops
awk '{print $2}' temp.cprops > temp.asm

bash ${pipeline}/finalize/construct-fasta-from-asm.sh temp.cprops temp.asm `basename $link .gz` | awk 'FILENAME==ARGV[1]{oname[FNR]=$1;next}$0~/>/{counter++;$0=">"oname[counter]}1' temp.cprops - > AaegL2.fasta

rm `basename $link .gz` `basename $link _genomic.fna.gz`"_assembly_report.txt" temp.cprops temp.asm