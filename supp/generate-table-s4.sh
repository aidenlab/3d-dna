#!/bin/bash
 
extract () {
	awk '{str[NR]=$0}END{print str[1]; print str[2]; print str[3]; print str[7]; print str[8]}' $1 | awk '{$1=$1}1'
}

in_chrom=`extract chr-length.fasta.stats | awk 'NR==1{print $NF}'`
		
	awk -v in_chrom=${in_chrom} '{tmp=$NF; $NF=""; print $0"\t"tmp}NR==1{tot_len=tmp}NR==5||NR==10||NR==15{pct=in_chrom/tot_len; print "In chromosome-length scaffolds\t"pct}' <(extract chr-length-and-small.fasta.stats) > S4.dat
	