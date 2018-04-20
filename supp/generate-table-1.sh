#!/bin/bash
 
extract () {
	awk '{str[NR]=$0}END{print str[1]; print str[2]; print str[3]; print str[7]; print str[8]}' $1 | awk '{$1=$1}1'
}
		
	awk '{tmp=$NF; $NF=""; print $0"\t"tmp}NR==5||NR==10||NR==15{print ""}' <(extract draft.fasta.stats && extract chr-length.fasta.stats && extract small.fasta.stats && extract tiny.fasta.stats) > T1.dat
	