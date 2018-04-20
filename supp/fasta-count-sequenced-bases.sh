#!/bin/bash

USAGE="./fasta-count-sequenced-bases.sh {-c <chrom_number>} <path_to_fasta>"

# handle options
while getopts "c:h" opt; do
case $opt in
	h) echo "$USAGE"
		exit 0
	;;
	c) re='^[-0-9]+$'
		if ! [[ $OPTARG =~ $re && $OPTARG -gt 0 ]] ; then
			echo ":( Error: Wrong syntax for number of chromosomes: counting all sequenced bases." >&2
		else
			echo ":) -c flag was triggered. Counting sequenced bases in the first $OPTARG scaffolds." >&1
			chrom=$OPTARG
		fi
	;;
	*) echo "$USAGE" >&1
		exit 1
	;;
esac
done

shift $(( OPTIND-1 ))


# handle arguments
fasta_file=$1

if [ -z $chrom ]
then
	parallel --pipepart -a ${fasta_file} --will-cite "awk '\$0!~/>/{gsub(\"N||n\",\"\"); c+=length}END{print c}'" | awk '{c+=$0}END{printf("Total seq bases: %'"'"'d\n", c)}'
else
	#echo "Individual chromosomes:"
	parallel --pipepart -a ${fasta_file} -k --will-cite "awk 'BEGIN{name=\"prev\"}\$0~/>/{if(c){print name, c}; name=\$1; c=0;next}{gsub(\"N||n\",\"\"); c+=length}END{print name, c}'" | awk -v chrom=${chrom} 'NR==1{name=$1; len=$2; next}$1=="prev"{len+=$2;next}{print name, len; name=$1; len=$2}END{print name, len}' | awk -v chrom=${chrom} 'NR<=chrom{c+=$2}{d+=$2}END{printf("Total chrom seq bases: %'"'"'d (%% of of total: %'"'"'f%%)\n", c, (c/d*100))}'
	
	#| tee >(head -n ${chrom}) '
fi
#printf("Total chrom seq bases: %'"'"'d, %% of total: %'"'"'d%%\n", c, (c/d*100))
