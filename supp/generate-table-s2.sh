#!/bin/bash
 
extract () {
	awk '{str[NR]=$0}END{print str[1]; print str[7]; print str[8]; print str[9]}' $1 | awk '{$1=$1}1'
}
		
	awk '{tmp=$NF; $NF=""; print $0"\t"tmp}NR==1{tot_seq=tmp}NR==5{print "% of Total Sequenced Base Pairs\t"tmp/tot_seq}NR==9{tot_attempt=tmp; print "% of Total Sequenced Base Pairs\t"tmp/tot_seq}NR==13{print "% of Total Sequenced Base Pairs\t"tmp/tot_seq}' <(extract input.fasta.stats && extract chr-length-and-small.fasta.stats && extract chr-length.fasta.stats && extract small.fasta.stats) > S2.dat
	