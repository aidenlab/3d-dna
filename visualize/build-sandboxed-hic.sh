#!/bin/bash

#### Description: Script to build final gapped assemblies (_HiC.assembly) in standard Juicebox format, with data 'sandboxed' to individual scaffolds and only chromosomes included as part of the hic file.
#### Usage: bash ./build-sandboxed-hic.sh [options] <path_to_HiC_assembly> <path_to_draft_mnd_file>.
#### Input: _HiC.assembly, draft mnd.txt file.
#### Output: hic file.
#### Parameters: number of chromosomes (default: all scaffolds in _HiC.assembly file), mapq threshold (0, 1, 30, all; default is 1 for mapq>=1).
#### Dependencies: Java; GNU Parallel.
#### Written by: Olga Dudchenko, version date 04/29/2020

USAGE="
*****************************************************
Visualizing draft genomes in juicebox: 18 July 2016

USAGE: ./build-sandboxed-hic.sh [options] <path_to_HiC_assembly> <path_to_draft_mnd_file>

DESCRIPTION:
This is a script to build final gapped assemblies (_HiC.assembly) in standard Juicebox format, with data 'sandboxed' to individual scaffolds and only chromosomes included as part of the hic file.

ARGUMENTS:
path_to_HiC_assembly
						Specify _HiC.assembly file (after adding all the gaps).
path_to_draft_mnd_file
						Specify path to mnd file describing pairwise Hi-C contacts between original assembly sequences.

OPTIONS:
-c chromosome_count
						Treat only the first n scaffolds in the _HiC.assembly as chromosomes.
-q mapq
						Build map for a specific mapq threshold (default is 1). 

**unprompted**
-r resolution_list
						Build for a specific list of resolutions (default is -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,1000).
-a path_to_draft_inter.txt
						Add metadata lifted from the original metadata file.
-m path_to_hic_mnd
						Use listed file and skip the remapping step.
-n
						Skip normalization.
-k
						Keep intermediate files.
-i
						Ignore mapq suffix.
-h
						Shows this help.
*****************************************************
"
## Defaults
mapq=1
use_parallel=true
res_string="2500000,1000000,500000,250000,100000,50000,25000,10000,5000,1000"

skip_norm=false
keep_files=false
verbose_mnd=0
ignore_mapq_suffix=false;

## HANDLE OPTIONS
while getopts "q:p:c:m:l:r:inkah" opt; do
case $opt in
	h) echo "$USAGE"
		exit 0
	;;
	c)  re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]] && [[ $OPTARG -gt 0 ]]; then
            echo ":) -c flag was triggered, treating the first ${OPTARG} scaffolds as chromosomes." >&1
            chrom_count=$OPTARG
        else
            echo ":( Wrong syntax for chromosome count. Exiting!" >&2
            exit 1
        fi
	;;
	q)  re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]]; then
            echo ":) -q flag was triggered, starting calculations for $OPTARG threshold mapping quality" >&1
            mapq=$OPTARG
        else
            echo ":( Wrong syntax for mapping quality. Exiting!" >&2
            exit 1
        fi
    ;;
	p)	if [ $OPTARG == true ] || [ $OPTARG == false ]; then
    	    echo ":) -p flag was triggered. Running with GNU Parallel support parameter set to $OPTARG." >&1
			use_parallel=$OPTARG
    	else
    		echo ":( Unrecognized value for -p flag. Exiting!" >&2
    		exit 1
    	fi
    ;;
	m) if [ -s $OPTARG ]; then 
			echo ":) Skipping remap step and using $OPTARG as premapped input." >&1
			remapped_mnd=$OPTARG
		else
			echo ":( Tentative remapped file not found. Exiting!" >&2
			exit 1
		fi		
	;;
	i)  ignore_mapq_suffix=true;
		echo ":) -i flag was triggered, building mapq without." >&1
	;;
	n)  skip_norm=true
		echo ":) -n flag was triggered, building maps without normalization." >&1
	;;
	k)  keep_files=true
		verbose_mnd=1
		echo ":) -k flag was triggered, will keep temporary files after completion." >&1
	;;
	r)  re='^[0-9]*(\,?[0-9]*)*$'
        if [[ $OPTARG =~ $re ]]; then
            echo ":) -r flag was triggered, starting calculations for resolution list: ${OPTARG}." >&1
            res_string=$OPTARG
        else
            echo ":( Wrong syntax for resolution flag. Exiting!" >&2
            exit 1
        fi
    ;;
    a)  if [ -f $OPTARG ]; then
    		echo ":) -a flag was triggered. Will add metadata based on file ${OPTARG}." >&1
    		original_metadata=$OPTARG
    	else
    		":( Tentative stats file not found. Exiting!" >&2
    		exit 1
		fi
	;;
	*)  echo ":( Illegal options. Exiting."
		echo "$USAGE"
		exit 1
	;;
esac
done

shift $(( OPTIND-1 ))

IFS=',' read -r -a res <<< "$res_string"

## HANDLE ARGUMENTS, TODO: check formats
assembly=$1
mnd=$2
[ $# -eq 2 ] && [ -s ${assembly} ] && [ -s ${mnd} ] || {	echo ":( Not sure how to parse your input or input files not found at intended locations. Exiting!" && echo "$USAGE" && exit 1 ; }

[ -z ${chrom_count} ] && echo >&2 ":| Warning: no explicit chromosome number listed, will attempt to build a hic map for all output scaffolds." && chrom_count=$(awk '$0!~/^>/{c++}END{print c}' ${assembly})

## CHECK DEPENDENCIES 

pipeline=`cd "$( dirname $0)" && cd .. && pwd`

type java >/dev/null 2>&1 || { echo >&2 ":( Java is not available, please install/add to path Java to run Juicer and Juicebox. Exiting!"; exit 1; }

if [ $use_parallel == true ]; then
	type parallel >/dev/null 2>&1 || { echo >&2 ":( GNU Parallel is not in the path. Continuing without parallelization."; use_parallel=false; }
fi


# ~~~~~~ ACTUAL SCRIPT ~~~~~~ #
genomeid=$(basename $assembly .assembly)

if [ -z ${remapped_mnd} ]; then
	cmd="awk -v label=\":::\" -v verbose=${verbose_mnd} -f ${pipeline}/edit/edit-mnd-according-to-new-assembly.awk ${assembly} - | awk -v sandbox=1 -v label=\"\" -v verbose=${verbose_mnd} -f ${pipeline}/lift/lift-input-mnd-to-assembly-mnd.awk ${assembly} -"

	if [ $use_parallel == false ]; then
		eval "cat $mnd | $cmd | sort --parallel=48 -S 32G -k2,2d -k6,6d | awk '{\$2=\"HiC_scaffold_\"\$2; \$6=\"HiC_scaffold_\"\$6}1'" > "temp."${genomeid}".mnd.txt"
	else
		#parallel -a $mnd --will-cite --jobs 80% --pipepart --block 1G "$cmd | sort -k2,2n -k6,6n" | sort -m -k2,2n -k6,6n --parallel=48 -S32G | awk '{$2="HiC_scaffold_"$2; $6="HiC_scaffold_"$6}1' > "temp."${genomeid}".mnd.txt"
		parallel -a $mnd --will-cite --jobs 80% --pipepart --block 1G "$cmd" | sort -k2,2n -k6,6n --parallel=48 -S32G | awk '{$2="HiC_scaffold_"$2; $6="HiC_scaffold_"$6}1' > "temp."${genomeid}".mnd.txt"
	fi
	remapped_mnd="temp."${genomeid}".mnd.txt"
fi

bash ${pipeline}/visualize/juicebox_tools.sh pre -q ${mapq} ${remapped_mnd} ${genomeid}${mapqsuf}.hic <(eval "awk -v n=${chrom_count} '\$0~/^>/{len[\$2]=\$3;len[-\$2]=\$3;next}{counter++}counter<=n{c=0; for(i=1;i<=NF;i++){c+=len[\$i]}; print \"HiC_scaffold_\"counter\"\t\"c}' ${assembly}")

[ $keep_files == true ] || rm "temp."${genomeid}".mnd.txt"
