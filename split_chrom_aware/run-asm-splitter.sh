#!/bin/bash

#### Description: Wrapper script to split megascaffold into individual chromosome-length scaffolds.
#### Written by: Sanjit Batra - sanjit.batra@bcm.edu. Version date 12/19/2016.

pipeline=`cd "$( dirname $0)" && cd .. && pwd`

USAGE="
*****************************************************
This is a wrapper script to split megascaffold into individual chromosome-length scaffolds. The output is a new asm file with additional chrom_num-1 number of lines where chrom_num is the number of chromosomes that are expected.

USAGE: ./run-asm-splitter.sh [options] <path_to_cprops> <path_to_asm> <path_to_contig_mnd>

DESCRIPTION:
This is a script to split an assembly chromosome into individual chromosomes.

ARGUMENTS:
path_to_cprops				Specify cprops file path.
path_to_asm					Specify asm file.
path_to_contig_mnd			Specify path to mnd file describing pairwise Hi-C contacts between assembly contigs.

OPTIONS:
-c chr_num					Number of chromosomes in the input genome (default is 23). 
-r true/false				Specify if the input is a rabl genome (true) or not (false), (default is false).
-h							Shows this help
*****************************************************
"

## Defaults
chr_num=23
rabl="false"

## unprompted
mapq=1

# non-rable default params
res=100000

# rable default params
rabl_wide=500000
rabl_narrow=50000

## HANDLE OPTIONS
while getopts "c:r:h" opt; do
case $opt in
	h) echo "$USAGE"
		exit 0
	;;
    c)  re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]] && [[ $OPTARG -gt 0 ]]; then
            echo ":) -n flag was triggered, running chromosome splitter with specified number of chromosomes = $OPTARG" >&1
            chr_num=$OPTARG
        else
            echo ":( Wrong syntax for specifying number of chromosomes. Using default number of chromosomes as ${chr_num}" >&2
        fi
	;;
	r)  if [ $OPTARG == true ] || [ $OPTARG == false ]; then
    	    echo ":) -r flag was triggered. Rabl status of genome is $OPTARG." >&1
			rabl=$OPTARG
    	else
    		echo ":( Unrecognized value for -r flag. Running with default parameters (-r ${rabl})." >&2
    	fi
    ;;
	*)  echo >&2 ":( Illegal options. Exiting."
		echo >&1 "$USAGE"
		exit 1
	;;
esac
done

shift $(( OPTIND-1 ))


## HANDLE ARGUMENTS, TODO: check formats
cprops=$1
asm=$2
mnd=$3

id=`basename ${cprops} .cprops`

[ -s ${cprops} ] && [ -s ${asm} ] && [ -s ${mnd} ] || {	echo ":( Not sure how to parse your input or input files not found at intended locations. Exiting!" && echo >&2 "$USAGE" && exit 1 ; }


## Handle zoom. TODO: get rid of zoom. What's the point if direct dump?

totlength=`awk '{total+=$3}END{print total}' ${cprops}`
scale=$(( 1 + $totlength / 2100000000 ))
if [ $scale -ne 1 ]; then
	res=$((res/scale))
	rabl_wide=$((rabl_wide/scale))
	rabl_narrow=$((rabl_narrow/scale))
fi

#Work with just the first line of the asm which is the first superscaffold
head -1 ${asm} > "asm_head"

maxl=`awk 'function abs(x){if(x<0){return -x}else{return x}}{if(FILENAME==ARGV[1]){clen[$2]=$3;next}}{if(FILENAME==ARGV[2]){for(i=1;i<=NF;i++){c+=clen[abs($i)]};next}}END{print c}' $cprops "asm_head"`
maxlength=$(( $maxl / $scale ))
buffer=10000000

#Lift over contig merged_nodups to assembly merged_nodups


bash ${pipeline}/lift/lift-input-mnd-to-asm-mnd.sh -s ${scale} ${cprops} "asm_head" ${mnd} > "asm_mnd.txt"

if [ $rabl == "false" ]; then

# Handle non-rabl situation

	#Dump matrix data
	bash ${pipeline}/split/mnd-dump.sh -q ${mapq} "asm_mnd.txt" ${res} > ${res}.jbdump

	#Run main recursive chrom splitter script
	python -W ignore ${pipeline}/split/recursive-chromosome-splitter.py ${res}.jbdump "boundaries.wig" ${chr_num} ${res}
	awk 'FNR>1{print $1}' "boundaries.wig" > "boundaries.list"

	#Here we perform two sanity checks on the boundaries: 
	#First, the number of lines in it should be chr_num : if not, we declare something went wrong during chromosome splitting
	nl=$(wc -l < "boundaries.wig")
	if [ $nl -ne $chr_num ]; then
		echo >&2 ":( Number of split chromosomes = $nl and input number of chromosomes = ${chr_num}, don't match. Refer to the hic map: continuing without splitting." && cp $asm ${id}.split.asm && exit
		#echo "Number of split chromosomes = "$nl
	fi
	#Second, we make sure that chromosome end positions are > 0
	test=`awk 'BEGIN{test=1}{if(\$1<=0){print ":( Chromosome boundary position not positive!" > "/dev/stderr"; test=0; exit}}END{print test}' "boundaries.wig"`
	[ $test -eq 0 ] && echo >&2 ":( Chromosome splitter failed. Refer to the hic map: continuing without splitting." && cp $asm ${id}.split.asm && exit
	
	#echo "Sanity check passed"

	#Produce boundary scaffolds
	awk -f ${pipeline}/split/find_scaffold_at_break.awk $cprops "asm_head" $mnd "boundaries.list" $scale | awk '{print $NF}' - > "scafs.list"	
	
	# cleanup
	rm ${res}.jbdump boundaries.wig boundaries.list
	
else

# Handle rabl situation

	#First we find the coarse chromosome boundaries using the coarse map juicebox dump
	bash ${pipeline}/split/mnd-dump.sh -q ${mapq} "asm_mnd.txt" ${rabl_wide} | awk '{if($1==0){print}}' - | sort -k 3,3nr | awk -v b=$buffer -v m=$maxlength '{if(!( ($2<b) || ($2>(m-b)) )){print $2}}' - | head -n$(( chr_num - 1 )) | sort -n > ${rabl_wide}.positions

	#Now we search for the maximum value in a radius of 1Mb??? 50kb centered at each of the coarse positions and find the scaffolds inside these positions lie	
	bash ${pipeline}/split/mnd-dump.sh -q ${mapq} "asm_mnd.txt" ${rabl_narrow} | awk '{if($1==0){print $2,$3}}' - | sort -k 1,1n | awk -f ${pipeline}/split/find_max.awk - ${rabl_wide}.positions | awk -f ${pipeline}/split/find_scaffold_at_rabl_peak.awk $cprops "asm_head" - $scale > "scafs.list"

	# cleanup

	rm ${rabl_wide}.positions

fi


# Finalize: break asm at scafs

awk -f ${pipeline}/split/split_asm_at_chr_boundaries.awk "scafs.list" "asm_head" > "chr.asm"
#cat "chr.asm" <(tail -n+2 $asm) > ${id}.split.asm
awk 'NR>1' $asm >> "chr.asm" && mv chr.asm ${id}.split.asm
awk '{$1=$1}1' ${id}.split.asm > ${id}.split.asm.tmp && mv ${id}.split.asm.tmp ${id}.split.asm

# cleanup
rm "scafs.list" "asm_head" "asm_mnd.txt"


