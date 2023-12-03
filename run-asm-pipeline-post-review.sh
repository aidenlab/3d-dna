#!/bin/bash
##########
#The MIT License (MIT)
#
# Copyright (c) 2018 Aiden Lab
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
##########
#
# 3D DNA de novo genome assembly pipeline: 180114 version.
#

echo `readlink -f $0`" "$*

USAGE_short="
*****************************************************
3D de novo assembly: version 180114

USAGE: ./run-asm-pipeline-post-review.sh [options] -r <review.assembly> <path_to_input_fasta> <path_to_input_mnd> 

DESCRIPTION:
This is a script to finalize assemblies (represented in input by draft fasta and deduplicated list of alignments of Hi-C reads to this fasta as produced by the Juicer pipeline) into chromosome-length fasta sequences, after review in Juicebox Assembly Tools Module (represented by review.assembly file). The script will produce an output fasta file, a Hi-C map of the final assembly, and a few supplementary annotation files to help review the result in Juicebox.

ARGUMENTS:
path_to_input_fasta			Specify file path to draft assembly fasta file.
path_to_input_mnd			Specify path to deduplicated list of alignments of Hi-C reads to the draft assembly fasta as produced by the Juicer pipeline: the merged_nodups file (mnd).

OPTIONS:
-r|--review path_to_review_assembly
			Path to review \".assembly\" file.
-i|--input input_size
			Specifies threshold input contig/scaffold size (default is 15000). Contigs/scaffolds smaller than input_size are going to be ignored. Only matters if running including seal. Should be the same as used for running the original script.
-s|--stage stage
			Assembly steps to run on top of the reviewed assembly, can be seal and finalize. Default is finalize.
-q|--mapq mapq
			Mapq threshold for final map(s) visualization, default is 1.
-g|--gap_size gap_size
			Gap size to be added between scaffolded sequences in the final chromosome-length scaffolds (default is 500).	
--sort-output
			Option to sort the chromosome-length scaffolds by size, in the descending order.
--build-gapped-map
			Option to output an additional contact map corresponding to the assembly after the gaps have been added between scaffolded sequences.
-h|--help			Shows this help. Type --help for a full set of options.

*****************************************************
"
pipeline=`cd "$( dirname $0)" && pwd`

## default parameter setup

input_size=15000 # contigs/scaffolds smaller than input_size are ignored
mapq=1	# default read mapping quality threshold for Hi-C scaffolder
gap_size=500	# default length of gaps to be added between scaffolded sequences in the chrom-length scaffolds
stage="finalize"	# by default run only final  pipeline
sort_output=false
build_gapped_map=false

############### HANDLE OPTIONS ###############

while :; do
	case $1 in
		-h|--help)
			echo "$USAGE_short" >&1
			exit 0
        ;;
        -r|--review) OPTARG=$2
			if [[ -f $OPTARG ]]; then
				echo " -r|--review flag was triggered, treating file $OPTARG as a JB4A review file for draft fasta in arguments." >&1
				review_assembly=$OPTARG
			else
				echo ":( File not found in the suggested review assembly file path. Exiting!" >&2
				exit 1
			fi
        	shift
        ;;
        -i|--input) OPTARG=$2
        	re='^[0-9]+$'
			if [[ $OPTARG =~ $re ]]; then
				echo " -i|--input flag was triggered, filtering draft contigs/scaffolds smaller than $OPTARG." >&1
				input_size=$OPTARG
			else
				echo ":( Wrong syntax for input size threshold. Using the default value ${input_size}." >&2
			fi
        	shift
        ;;
        -s|--stage) OPTARG=$2
			if [ "$OPTARG" == "seal" ] || [ "$OPTARG" == "finalize" ]; then
			echo " -s|--stage flag was triggered, fast-forwarding to \"$OPTARG\" pipeline section." >&1
			stage=$OPTARG
			else
				echo ":( Wrong syntax for pipeline stage. Exiting!" >&2
			fi
        	shift
        ;;
        -q|--mapq) OPTARG=$2
			re='^[0-9]+$'
			if [[ $OPTARG =~ $re ]]; then
				echo " -q|--mapq flag was triggered, scaffolding using reads with at least $OPTARG mapping quality." >&1
				mapq=$OPTARG
			else
				echo ":( Wrong syntax for mapping quality. Using the default value ${mapq}." >&2
			fi
        	shift
        ;;
        -g|--gap-size) OPTARG=$2
			re='^[0-9]+$'
			if [[ $OPTARG =~ $re ]]; then
				echo " -g|--gap-size flag was triggered, will add gaps of size $OPTARG between scaffolded sequences in the chromosome-length scaffolds." >&1
				gap_size=$OPTARG
			else
				echo ":( Wrong syntax for gap size parameter value. Using the default value ${gap_size}." >&2
			fi
        	shift 
		;;
## long menu options
		--sort-output)
			echo " --sort-output was triggered, will sort output scaffolds by size." >&1
			sort_output=true
		;;
		--build-gapped-map)
			echo " --build-gapped-map was triggered, will build an additional hic file corresponding to final assembly with gaps added between draft sequences." >&1
			build_gapped_map=true
		;;
		--) # End of all options
			shift
			break
		;;
		-?*)
			echo ":| WARNING: Unknown option. Ignoring: ${1}" >&2
		;;
		*) # Default case: If no more options then break out of the loop.
			break
	esac
	shift
done

############### HANDLE EXTERNAL DEPENDENCIES ###############

##	GNU Parallel Dependency
parallel="false"
if hash parallel 2>/dev/null; then
        ver=`parallel --version | awk 'NR==1{print \$3}'`
        [ $ver -ge 20150322 ] && parallel="true"
fi
[ $parallel == "false" ] && echo ":| WARNING: GNU Parallel version 20150322 or later not installed. We highly recommend to install it to increase performance. Starting pipeline without parallelization!" >&2

############### HANDLE ARGUMENTS ###############

[ -z ${review_assembly} ] || [ -z $1 ] || [ -z $2 ] && echo >&2 "Not sure how to parse your input: files not listed or not found at expected locations. Exiting!" && echo >&2 "$USAGE_short" && exit 1

[ ! -s ${review_assembly} ] || [ ! -s $1 ] || [ ! -s $2 ] && echo >&2 "Not sure how to parse your input: files not listed or not found at expected locations. Exiting!" && echo >&2 "$USAGE_short" && exit 1

## TODO: check file format

if [ "$#" -ne 2 ]; then
    echo >&2 "Illegal number of arguments. Please double check your input. Exiting!" && echo >&2 "$USAGE_short" && exit 1
fi

orig_fasta=$1
orig_mnd=$2

genomeid=$(basename "$orig_fasta" .fasta)
genomeid=$(basename "$genomeid" .fna)
genomeid=$(basename "$genomeid" .fa)

if [ ${orig_mnd} != ${genomeid}.mnd.txt ]; then
	if [ -f ${genomeid}.mnd.txt ]; then
		cmp --silent ${orig_mnd} ${genomeid}.mnd.txt || { echo >&2 ":( Please remove or rename file ${genomeid}.mnd.txt. Exiting!" && exit 1; }
	fi
	ln -sf ${orig_mnd} ${genomeid}.mnd.txt
fi

orig_mnd=${genomeid}.mnd.txt

## TODO: add some checks to make more user-proof against swapping
[ -f ${genomeid}.cprops ] || awk -f ${pipeline}/utils/generate-sorted-cprops-file.awk ${orig_fasta} > ${genomeid}.cprops

orig_cprops=${genomeid}.cprops
[ ! -f ${orig_cprops} ] && echo >&2 ":( No cprops file found. Please rerun the pipeline from scratch. Exiting!" && exit 1

##	calculate zoom
# TODO: move this to mismatch detector, pass only scale
totlength=`awk '{total+=$3}END{print total}' ${orig_cprops}`
scale=$(( 1 + $totlength / 2100000000 ))
if [ $scale -ne 1 ]; then
	editor_coarse_resolution=$((editor_coarse_resolution/scale))
	editor_coarse_region=$((editor_coarse_region/scale))
	editor_fine_resolution=$((editor_fine_resolution/scale))
	
	polisher_coarse_resolution=$((polisher_coarse_resolution/scale))
	polisher_coarse_region=$((polisher_coarse_region/scale))
	polisher_fine_resolution=$((polisher_fine_resolution/scale))
	
	splitter_coarse_resolution=$((splitter_coarse_resolution/scale))
	splitter_coarse_region=$((splitter_coarse_region/scale))
	splitter_fine_resolution=$((splitter_fine_resolution/scale))
fi


if [ "$stage" != "finalize" ]; then

############### SEALING ###############


	awk -v cprops=${genomeid}.split.cprops -v asm=${genomeid}.split.asm '$1~/^>/{$1=substr($1,2); print > cprops;next}{print > asm}' ${review_assembly}	
	add_options=""
	[ "$sort_output" == "true" ] && add_options="--sort-output"
	[ "$build_gapped_map" == "true" ] && add_tions=${add_options}" --build-gapped-map"
	if [ "$add_options" != "" ]; then
		bash ${pipeline}/run-asm-pipeline.sh -s seal -i ${input_size} -g ${gap_size} ${add_options} ${orig_fasta} ${orig_mnd}
	else
		bash ${pipeline}/run-asm-pipeline.sh -s seal -i ${input_size} -g ${gap_size} ${orig_fasta} ${orig_mnd}
	fi
	
else

############### FINALIZING ###############

	echo "###############" >&1
	echo "Finilizing output:" >&1
	
	awk -v cprops=${genomeid}.final.cprops -v asm=${genomeid}.final.asm '$1~/^>/{$1=substr($1,2); print > cprops;next}{print > asm}' ${review_assembly}	

	if [ "$sort_output" == "true" ]; then
		awk -v input_size=${input_size} -f ${pipeline}/utils/sort-asm-by-size.awk ${genomeid}.final.cprops ${genomeid}.final.asm > ${genomeid}.final.asm.tmp && mv ${genomeid}.final.asm.tmp ${genomeid}.final.asm
	fi

	# build final map
	bash ${pipeline}/edit/edit-mnd-according-to-new-cprops.sh ${genomeid}.final.cprops ${orig_mnd} > ${genomeid}.final.mnd.txt
	bash ${pipeline}/visualize/run-asm-visualizer.sh -p ${parallel} -q ${mapq} -i -c ${genomeid}.final.cprops ${genomeid}.final.asm ${genomeid}.final.mnd.txt
	rm ${genomeid}.final.mnd.txt

	# build final fasta
	#awk -f ${pipeline}/edit/edit-fasta-according-to-new-cprops.awk ${genomeid}.final.cprops ${orig_fasta} > ${genomeid}.final.fasta
    python3 ${pipeline}/edit/edit-fasta-according-to-new-cprops.py ${genomeid}.final.cprops ${orig_fasta} | seqkit seq -w 80 > ${genomeid}.final.fasta
 
	bash ${pipeline}/finalize/finalize-output.sh -s ${input_size} -l ${genomeid} -g ${gap_size} ${genomeid}.final.cprops ${genomeid}.final.asm ${genomeid}.final.fasta final

	# if requested build HiC map with added gaps
	if [ "$build_gapped_map" == "true" ]; then
		awk -f ${pipeline}/utils/convert-assembly-to-cprops-and-asm.awk ${genomeid}.FINAL.assembly
		bash ${pipeline}/edit/edit-mnd-according-to-new-cprops.sh ${genomeid}.FINAL.cprops ${orig_mnd} > ${genomeid}.FINAL.mnd.txt
		bash ${pipeline}/visualize/run-assembly-visualizer.sh -p ${parallel} -q ${mapq} -i -c ${genomeid}.FINAL.assembly ${genomeid}.FINAL.mnd.txt
		rm ${genomeid}.FINAL.mnd.txt ${genomeid}.FINAL.cprops ${genomeid}.FINAL.asm
	fi
fi
