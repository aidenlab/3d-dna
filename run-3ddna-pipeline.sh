#!/bin/bash
##########
#The MIT License (MIT)
#
# Copyright (c) 2021 Aiden Lab
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
# 3D-DNA de novo genome assembly pipeline.
#

version=210623

echo `readlink -f $0`" "$*

echo "version: "${version}

#set -x

USAGE_short="
*****************************************************
3D de novo assembly: version "$version"

USAGE: ./run-3ddna-pipeline.sh [options] <path_to_input_fasta> <path_to_input_mnd> 

DESCRIPTION:
This is a script to assemble draft assemblies (represented in input by draft fasta and deduplicated list of alignments of Hi-C reads to this fasta as produced by the Juicer pipeline) into chromosome-length scaffolds. The script will produce an output fasta file, a Hi-C map of the final assembly, and a few supplementary annotation files to help review the result in Juicebox.

ARGUMENTS:
path_to_input_fasta
							Specify file path to draft assembly fasta file.
path_to_input_mnd
							Specify path to deduplicated list of alignments of Hi-C reads to the draft assembly fasta as produced by the Juicer pipeline: the merged_nodups file (mnd).

OPTIONS:
-r|--review path_to_review_assembly_file
							Specifies a path to a JBAT review.assembly file.
-c|--chrom_count number_of_chromosomes
							Instructs the pipeline to treat the first <number_of_chromosomes> Hi-C scaffolds as chromosomes. Usually used in conjunction with -r once the scaffolds have been manually reviewed and polished. The flag will trigger generating a standard sandboxed hic contact file fort the chromosome-length assembly.
--organism organism_name
							Add an [organism=\"organism_name\"] tag to output sequence names.
--isolate isolate_name
							Add a [isolate=\"isolate_name\"] tag to output sequence names.
-h|--help
							Type -h to shows this help. Type --help for a full set of options.
*****************************************************
"

USAGE_long="
*****************************************************
3D de novo assembly: version "$version"

USAGE: ./run-3ddna-pipeline.sh [options] <path_to_input_fasta> <path_to_input_mnd> 

DESCRIPTION:
This is a script to assemble draft assemblies (represented in input by draft fasta and deduplicated list of alignments of Hi-C reads to this fasta as produced by the Juicer pipeline) into chromosome-length scaffolds. The script will produce an output fasta file, a Hi-C map of the final assembly, and a few supplementary annotation files to help review the result in Juicebox.

ARGUMENTS:
path_to_input_fasta
			Specify file path to draft assembly fasta file.
path_to_input_mnd
			Specify path to deduplicated list of alignments of Hi-C reads to the draft assembly fasta as produced by the Juicer pipeline: the merged_nodups file (mnd).

OPTIONS:
-h
			Shows main options.
--help
			Shows this help.

-r|--review path_to_review_assembly_file
			Specify a path to a JBAT review.assembly file.
-c|--chromosomes number_of_chromosomes
			Instructs the pipeline to treat the first <number_of_chromosomes> Hi-C scaffolds as chromosomes. Usually used in conjunction with -r once the scaffolds have been manually reviewed and polished. The flag will trigger generating a standard sandboxed hic contact file fort the chromosome-length assembly. TODO: add appropriate labels to fasta.

-i|--input input_size
			Specifies threshold input contig/scaffold size (default is 15000). Contigs/scaffolds smaller than input_size are going to be ignored.
-e|--edits number_of_edit_rounds
			Specifies number of iterative rounds for misjoin correction (default is 2).
-s|--stage stage
			Fast forward to later assembly steps, can be polish, split, seal, merge and finalize.
--organism organism_name
			Add an [organism=\"organism_name\"] tag to output sequence names.
--isolate isolate_name
			Add a [isolate=\"isolate_name\"] tag to output sequence names.

ADDITIONAL OPTIONS:
**scaffolder**
-q|--mapq mapq					Mapq threshold for scaffolding and visualization (default is 1).

**misjoin detector**
--editor-coarse-resolution editor_coarse_resolution
			Misjoin editor coarse matrix resolution, should be one of the following: 2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000, 1000 (default is 25000).
--editor-coarse-region editor_coarse_region
			Misjoin editor triangular motif region size (default is 125000).
--editor-coarse-stringency editor_coarse_stringency
			Misjoin editor stringency parameter (default is 55).
--editor-saturation-centile editor_saturation_centile
			Misjoin editor saturation parameter (default is 5).
--editor-fine-resolution editor_fine_resiolution
			Misjoin editor fine matrix resolution, should be one of the following: 2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000, 1000 (default is 1000).
--editor-repeat-coverage editor_repeat_coverage
			Misjoin editor threshold repeat coverage (default is 2). 

**polisher**
--polisher-input-size polisher_input_size
			Polisher input size threshold. Scaffolds smaller than polisher_input_size are going to be placed into unresolved (default is 1000000).
--polisher-coarse-resolution editor_coarse_resolution
			Polisher coarse matrix resolution, should be one of the following: 2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000, 1000 (default is 25000).
--polisher-coarse-region editor_coarse_region
			Polisher  triangular motif region size (default is 3000000).
--polisher-coarse-stringency editor_coarse_stringency
			Polisher stringency parameter (default is 55).
--polisher-saturation-centile editor_saturation_centile
			Polisher saturation parameter (default is 5).
--polisher-fine-resolution editor_fine_resiolution
			Polisher fine matrix resolution, should be one of the following: 2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000, 1000 (default is 1000).

**splitter**
--splitter-input-size splitter_input_size
			Splitter input size threshold. Scaffolds smaller than polisher_input_size are going to be placed into unresolved (Default: 1000000).
--splitter-coarse-resolution splitter_coarse_resolution
			Splitter coarse matrix resolution, should be one of the following: 2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000, 1000 (Default: 25000).
--splitter-coarse-region splitter_coarse_region
			Splitter  triangular motif region size (Default: 3000000).
--splitter-coarse-stringency splitter_coarse_stringency
			Splitter stringency parameter (Default: 55).
--splitter-saturation-centile splitter_saturation_centile
			Splitter saturation parameter (Default: 5).
--splitter-fine-resolution splitter_fine_resiolution
			Splitter fine matrix resolution, should be one of the following: 2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000, 1000 (Default: 1000).

**merger**
--merger-search-band merger_search_band		
			Distance (in bp) within which to locally search for alternative haplotypes to a given contig or scaffold, from the position of their suggested incorporation in the assembly. The larger the original input contigs/scaffolds, the larger band size it might be necessary to set. Default: 3000000.
--merger-alignment-score merger_alignment_score
			Minimal LASTZ alignment score for nearby sequences (located in the assembly within the distance defined by the merger_search_band parameter) to be recongnized as alternative haplotypes. Default: 50000000.
--merger-alignment-identity merger_alignment_identity
			Minimal identity score required from similar nearby sequences (per length) for them to be classified as alternative haplotypes. Default: 20.
--merger-alignment-length	merger_alignment_length
			Minimal length necessary to recognize similar nearby sequences as alternative haplotypes. Default: 20000.
--merger-lastz-options	merger_lastz_options
			Option string to customize LASTZ alignment. Default: \"--gfextend\ --gapped\ --chain=200,200\"		

**finalizer**
-g|--gap-size gap_size
			Gap size to be added between scaffolded sequences in the final chromosome-length scaffolds (default is 500).

**supplementary**
--early-exit
			Option to exit after first round of scaffolding.
--fast-start
			Option to pick up processing assuming the first round of scaffolding is done. In conjunction with --early-exit this option is to help tune the parameters for best performance.
--sort-chromosomes [number_of_chromosomes_to_sort]
			Option to sort the chromosome-length scaffolds by size, in the descending order. If a value is provided with the flag only the first <number_of_chromosomes_to_sort> will be sorted by size. Otherwise all scaffolds will be sorted.
*****************************************************
"

pipeline=`cd "$( dirname $0)" && pwd`

## default parameter setup

input_size=15000 # contigs/scaffolds smaller than input_size are ignored
MAX_ROUNDS=2	# usually 2 is adequate

mapq=1	# default read mapping quality threshold for Hi-C scaffolder

# misassembly detector and editor default params
editor_saturation_centile=5
editor_coarse_resolution=25000
editor_coarse_region=125000
editor_coarse_stringency=55
editor_coarse_norm="KR"	
editor_fine_resolution=1000
editor_repeat_coverage=2

# polisher default params
polisher_coarse_resolution=100000	
polisher_fine_resolution=1000
polisher_coarse_region=3000000
polisher_coarse_stringency=55
polisher_saturation_centile=5
polisher_input_size=1000000

# splitter detection default params
splitter_input_size=100000 # in principle don't really need this, just moves smallish scaffolds to the back
splitter_coarse_resolution=100000	
splitter_fine_resolution=1000
splitter_coarse_region=3000000
splitter_coarse_stringency=55
splitter_saturation_centile=5

# finalizer default params
gap_size=500	# default length of gaps to be added between scaffolded sequences in the chrom-length scaffolds

# supplementary options
stage=""	# by default run full pipeline
early=false
fast=false
sort_chromosomes=1 # by default do not sort chromosomes

############### HANDLE OPTIONS ###############

while :; do
	case $1 in
		-h)
			echo "$USAGE_short" >&1
			exit 0
        ;;
  		--help)
			echo "$USAGE_long" >&1
			exit 0
        ;;

## SHORT MENU OPTIONS
        -r|--review) OPTARG=$2
			if [[ -s $OPTARG ]]; then
				echo " -r|--review flag was triggered, treating file $OPTARG as a JBAT review file for draft fasta in arguments." >&1
				review_assembly=$OPTARG
			else
				echo ":( File not found in the suggested review assembly file path or file is empty. Exiting!" >&2
				exit 1
			fi
        	shift
        ;;        
		-c|--chromosomes) OPTARG=$2
			re='^[0-9]+$'
			if [[ $OPTARG =~ $re ]]; then
					echo " -c|--chromosomes flag was triggered, will build an additional standard \"sandboxed\" hic file corresponding to _HiC.assembly with $OPTARG chromosome(s)." >&1
					chromosome_count=$OPTARG
			else
					echo ":( Wrong syntax for chromosome count parameter value. Exiting!" >&2
					exit 1
			fi
			shift
        ;;
        -e|--edits) OPTARG=$2
        	re='^[0-9]+$'
			if [[ $OPTARG =~ $re ]] && [[ $OPTARG -ge 0 ]]; then
				echo " -r|--rounds flag was triggered, will run $OPTARG round(s) of misjoin correction." >&1
				MAX_ROUNDS=$OPTARG
			else
				echo ":( Wrong syntax for number of iterative rounds of misjoin correction. Using the default value ${MAX_ROUNDS}." >&2
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
# scaffolder
        -q|--mapq) OPTARG=$2 ##TODO: check that propagates consistently, not tested sufficiently
			re='^[0-9]+$'
			if [[ $OPTARG =~ $re ]]; then
				echo " -q|--mapq flag was triggered, scaffolding using reads with at least $OPTARG mapping quality." >&1
				mapq=$OPTARG
			else
				echo ":( Wrong syntax for mapping quality. Using the default value ${mapq}." >&2
			fi
        	shift
        ;;
# organizational
		-s|--stage) OPTARG=$2
			if [ "$OPTARG" == "scaffold" ] || [ "$OPTARG" == "polish" ] || [ "$OPTARG" == "split" ] || [ "$OPTARG" == "seal" ] || [ "$OPTARG" == "finalize" ]; then
				echo " -s|--stage flag was triggered, fast-forwarding to \"$OPTARG\" pipeline section." >&1
				stage=$OPTARG
			else
				echo ":( Wrong syntax for pipeline stage. Exiting!" >&2
			fi
        	shift
        ;;
        
## LONG MENU OPTIONS
# misjoin editor
        --editor-saturation-centile) OPTARG=$2
			re='^[0-9]+\.?[0-9]*$'
			if [[ $OPTARG =~ $re ]] && [[ ${OPTARG%.*} -ge 0 ]] && ! [[ "$OPTARG" =~ ^0*(\.)?0*$ ]] && [[ $((${OPTARG%.*} + 1)) -le 100 ]]; then
				echo " --editor-saturation-centile flag was triggered, misjoin editor saturation parameter set to ${OPTARG}%." >&1
				editor_saturation_centile=$OPTARG
			else
				echo ":( Wrong syntax for misjoin editor saturation threshold. Using the default value ${editor_saturation_centile}%." >&2
			fi
        	shift
        ;;
        --editor-coarse-resolution) OPTARG=$2
        	re='^[0-9]+$'	## TODO: specify/generalize re matrix resolutions size
			if [[ $OPTARG =~ $re ]]; then
				echo " --editor-coarse-resolution flag was triggered, misjoin editor coarse matrix resolution set to $OPTARG." >&1
				editor_coarse_resolution=$OPTARG
			else
				echo ":( Wrong syntax for misjoin editor coarse matrix resolution. Using the default value ${editor_coarse_resolution}." >&2
			fi
        	shift
        ;;
        --editor-coarse-region) OPTARG=$2
			re='^[0-9]+$'
			if [[ $OPTARG =~ $re ]]; then
				echo " --editor-coarse-region flag was triggered, misjoin editor coarse resolution depletion region size set to $OPTARG." >&1
				editor_coarse_region=$OPTARG
			else
				echo ":( Wrong syntax for misjoin editor coarse resolution depletion region size. Using the default value ${editor_coarse_region}." >&2
			fi
        	shift
        ;;
        --editor-coarse-stringency) OPTARG=$2
			re='^[0-9]+$'
			if [[ $OPTARG =~ $re ]] && [[ $OPTARG -gt 0 ]] && [[ $OPTARG -lt 100 ]]; then
				echo " --editor-coarse-stringency flag was triggered, misjoin detection stringency parameter set to $OPTARG%." >&1
				editor_coarse_stringency=$OPTARG
			else
				echo ":( Wrong syntax for misjoin detection stringency parameter. Using the default value ${editor_coarse_stringency}%." >&2
			fi
        	shift
        ;;
        --editor-coarse-norm) OPTARG=$2
        	if [ $OPTARG == NONE ] || [ $OPTARG == VC ] || [ $OPTARG == VC_SQRT ] || [ $OPTARG == KR ]; then
    	    	echo " --editor-coarse-norm flag was triggered. Type of norm chosen for the contact matrix is $OPTARG." >&1
				editor_coarse_norm=$OPTARG
    		else
    			echo ":( Unrecognized value for -b flag. Running with default parameters (-b NONE)." >&2
    		fi
        	shift
        ;;
        --editor-fine-resolution) OPTARG=$2
			re='^[0-9]+$'	## TODO: specify/generalize re matrix resolutions size
			if [[ $OPTARG =~ $re ]]; then
				echo " --editor-fine-resolution flag was triggered, misjoin detection fine matrix resolution set to $OPTARG." >&1
				editor_fine_resolution=$OPTARG
			else
				echo ":( Wrong syntax for misjoin editor fine matrix resolution. Using the default value ${editor_fine_resolution}." >&2
			fi
        	shift
        ;;
        --editor-repeat-coverage) OPTARG=$2
        	re='^[0-9]+$'
			if [[ $OPTARG =~ $re ]] && [[ $OPTARG -gt 0 ]]; then
				echo " --editor-repeat-coverage flag was triggered, threshold repeat coverage parameter set to $OPTARG." >&1
				editor_repeat_coverage=$OPTARG
			else
				echo ":( Wrong syntax for misjoin detection stringency parameter. Using the default value ${editor_repeat_coverage}." >&2
			fi
        	shift
        ;;
# polisher
        --polisher-input-size) OPTARG=$2
        	re='^[0-9]+$'
			if [[ $OPTARG =~ $re ]]; then
				echo " --polisher-input-size flag was triggered, excluding scaffolds smaller than $OPTARG when polishing." >&1
				polisher_input_size=$OPTARG
			else
				echo ":( Wrong syntax for polisher scaffold input size threshold. Using the default value ${polisher_input_size}." >&2
			fi
        	shift
        ;;
        --polisher-saturation-centile) OPTARG=$2
        	re='^[0-9]+\.?[0-9]*$'
			if [[ $OPTARG =~ $re ]] && [[ ${OPTARG%.*} -ge 0 ]] && ! [[ "$OPTARG" =~ ^0*(\.)?0*$ ]] && [[ $((${OPTARG%.*} + 1)) -le 100 ]]; then
				echo " --polisher-saturation-centile flag was triggered, polisher saturation parameter set to ${OPTARG}%." >&1
				polisher_saturation_centile=$OPTARG
			else
				echo ":( Wrong syntax for polisher saturation threshold. Using the default value ${polisher_saturation_centile}%." >&2
			fi
        	shift
        ;;
        --polisher-coarse-resolution) OPTARG=$2
        	re='^[0-9]+$'	## TODO: specify/generalize re matrix resolutions size
			if [[ $OPTARG =~ $re ]]; then
				echo " --polisher-coarse-resolution flag was triggered, polisher coarse matrix resolution set to $OPTARG." >&1
				polisher_coarse_resolution=$OPTARG
			else
				echo ":( Wrong syntax for polisher coarse matrix resolution. Using the default value ${polisher_coarse_resolution}." >&2
			fi
        	shift
        ;;
        --polisher-coarse-region) OPTARG=$2
			re='^[0-9]+$'
			if [[ $OPTARG =~ $re ]]; then
				echo " --polisher-coarse-region flag was triggered, polisher coarse resolution depletion region size set to $OPTARG." >&1
				polisher_coarse_region=$OPTARG
			else
				echo ":( Wrong syntax for polisher coarse resolution depletion region size. Using the default value ${polisher_coarse_region}." >&2
			fi
        	shift
        ;;
        --polisher-coarse-stringency) OPTARG=$2
			re='^[0-9]+$'
			if [[ $OPTARG =~ $re ]] && [[ $OPTARG -gt 0 ]] && [[ $OPTARG -lt 100 ]]; then
				echo " --polisher-coarse-stringency flag was triggered, polisher stringency parameter set to $OPTARG%." >&1
				polisher_coarse_stringency=$OPTARG
			else
				echo ":( Wrong syntax for polisher stringency parameter. Using the default value ${polisher_coarse_stringency}%." >&2
			fi
        	shift
        ;;
        --polisher-fine-resolution) OPTARG=$2
			re='^[0-9]+$'	## TODO: specify/generalize re matrix resolutions size
			if [[ $OPTARG =~ $re ]]; then
				echo " --polisher-fine-resolution flag was triggered, polisher fine matrix resolution set to $OPTARG." >&1
				polisher_fine_resolution=$OPTARG
			else
				echo ":( Wrong syntax for polisher fine matrix resolution. Using the default value ${polisher_fine_resolution}." >&2
			fi
        	shift
        ;;

# splitter
        --splitter-input-size) OPTARG=$2 ##TODO: should get rid of this, don't think I really need it
        	re='^[0-9]+$'
			if [[ $OPTARG =~ $re ]]; then
				echo " --splitter-input-size flag was triggered, excluding scaffolds smaller than $OPTARG when splitting." >&1
				splitter_input_size=$OPTARG
			else
				echo ":( Wrong syntax for splitter scaffold input size threshold. Using the default value ${splitter_input_size}." >&2
			fi
        	shift
        ;;
        --splitter-saturation-centile) OPTARG=$2
        	re='^[0-9]+\.?[0-9]*$'
			if [[ $OPTARG =~ $re ]] && [[ ${OPTARG%.*} -ge 0 ]] && ! [[ "$OPTARG" =~ ^0*(\.)?0*$ ]] && [[ $((${OPTARG%.*} + 1)) -le 100 ]]; then
				echo " --splitter-saturation-centile flag was triggered, splitter saturation parameter set to ${OPTARG}%." >&1
				splitter_saturation_centile=$OPTARG
			else
				echo ":( Wrong syntax for splitter saturation threshold. Using the default value ${splitter_saturation_centile}%." >&2
			fi
        	shift
        ;;
        --splitter-coarse-resolution) OPTARG=$2
        	re='^[0-9]+$'	## TODO: specify/generalize re matrix resolutions size
			if [[ $OPTARG =~ $re ]]; then
				echo " --splitter-coarse-resolution flag was triggered, splitter coarse matrix resolution set to $OPTARG." >&1
				splitter_coarse_resolution=$OPTARG
			else
				echo ":( Wrong syntax for splitter coarse matrix resolution. Using the default value ${splitter_coarse_resolution}." >&2
			fi
        	shift
        ;;
        --splitter-coarse-region) OPTARG=$2
			re='^[0-9]+$'
			if [[ $OPTARG =~ $re ]]; then
				echo " --splitter-coarse-region flag was triggered, splitter coarse resolution depletion region size set to $OPTARG." >&1
				splitter_coarse_region=$OPTARG
			else
				echo ":( Wrong syntax for splitter coarse resolution depletion region size. Using the default value ${splitter_coarse_region}." >&2
			fi
        	shift
        ;;
        --splitter-coarse-stringency) OPTARG=$2
			re='^[0-9]+$'
			if [[ $OPTARG =~ $re ]] && [[ $OPTARG -gt 0 ]] && [[ $OPTARG -lt 100 ]]; then
				echo " --splitter-coarse-stringency flag was triggered, splitter stringency parameter set to $OPTARG%." >&1
				splitter_coarse_stringency=$OPTARG
			else
				echo ":( Wrong syntax for splitter stringency parameter. Using the default value ${splitter_coarse_stringency}%." >&2
			fi
        	shift
        ;;
        --splitter-fine-resolution) OPTARG=$2
			re='^[0-9]+$'	## TODO: specify/generalize re matrix resolutions size
			if [[ $OPTARG =~ $re ]]; then
				echo " --splitter-fine-resolution flag was triggered, splitter fine matrix resolution set to $OPTARG." >&1
				splitter_fine_resolution=$OPTARG
			else
				echo ":( Wrong syntax for splitter fine matrix resolution. Using the default value ${splitter_fine_resolution}." >&2
			fi
        	shift
        ;;
# finalizer        
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
# supplementary options:
		--early-exit)
			echo " --early-exit flag was triggered, will do early exit." >&1
			early=true
		;;
		--fast-start)
			echo " --fast-start flag was triggered, will start assuming first iterative round and map are available." >&1
			fast=true
		;;
		--sort-chromosomes) OPTARG=$2
			re='^[0-9]+$'
			if [[ $OPTARG =~ $re ]]; then
				echo " --sort-chromosomes was triggered with $OPTARG parameter, will sort the first $OPTARG output scaffolds by size." >&1
				sort_chromosomes=$OPTARG
				shift
			else
				echo " --sort-chromosomes was triggered without an argument, will sort all output scaffolds by size." >&1
				sort_chromosomes=0
			fi
		;;
		--organism) OPTARG=$2
			echo " --organism flag was triggered with $OPTARG parameter, will add an [organism=$OPTARG] tag to output sequence names." >&1
			organism=$OPTARG
			shift
		;;
		--isolate) OPTARG=$2
			echo " --isolate flag was triggered with $OPTARG parameter, will add an [isolate=$OPTARG] tag to output sequence names." >&1
			isolate=$OPTARG
			shift
		;;
# TODO: merger, sealer, etc options              
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

## check parameters for compatibility
[[ ${editor_coarse_region} -le ${editor_coarse_resolution} ]] && echo >&2 ":( Requested depletion region size ${editor_coarse_region} and bin size ${editor_coarse_resolution} parameters for editor are incompatible. Run ${pipeline}/edit/run-mismatch-detector.sh -h for instructions. Exiting!" && exit 1
[[ ${editor_coarse_resolution} -le ${editor_fine_resolution} ]] && echo >&2 ":( Requested mismatch localization resolution ${editor_fine_resolution} and coarse search bin size ${editor_coarse_resolution} parameters for editor are incompatible. Run ${pipeline}/edit/run-mismatch-detector.sh -h for instructions. Exiting!" && exit 1

[[ ${polisher_coarse_region} -le ${polisher_coarse_resolution} ]] && echo >&2 ":( Requested depletion region size ${polisher_coarse_region} and bin size ${polisher_coarse_resolution} parameters for polisher are incompatible. Run ${pipeline}/edit/run-mismatch-detector.sh -h for instructions. Exiting!" && exit 1
[[ ${polisher_coarse_resolution} -le ${polisher_fine_resolution} ]] && echo >&2 ":( Requested mismatch localization resolution ${polisher_fine_resolution} and coarse search bin size ${polisher_coarse_resolution} parameters for polisher are incompatible. Run ${pipeline}/edit/run-mismatch-detector.sh -h for instructions. Exiting!" && exit 1

[[ ${splitter_coarse_region} -le ${splitter_coarse_resolution} ]] && echo >&2 ":( Requested depletion region size ${splitter_coarse_region} and bin size ${splitter_coarse_resolution} parameters for splitter are incompatible. Run ${pipeline}/edit/run-mismatch-detector.sh -h for instructions. Exiting!" && exit 1
[[ ${splitter_coarse_resolution} -le ${splitter_fine_resolution} ]] && echo >&2 ":( Requested mismatch localization resolution ${splitter_fine_resolution} and coarse search bin size ${splitter_coarse_resolution} parameters for splitter are incompatible. Run ${pipeline}/edit/run-mismatch-detector.sh -h for instructions. Exiting!" && exit 1

([[ ! -z ${review_assembly} ]] && ! ([ "$stage" != "" ] || [ "$stage" != "seal" ] || [ "$stage" != "finalize" ])) && echo >&2 ":( Review mode is incompatible with the set stage parameter $stage. Exiting!" && exit 1

([[ ! -z ${review_assembly} ]] && [[ "$stage" != "seal" ]]) && stage="finalize" # in review mode the stage can be either seal or finalize

([[ ${chromosome_count} -gt 1 ]] && [[ -z ${review_assembly} ]]) && echo >&2 ":| WARNING: -c|--chromosomes flag has been triggered while not in review mode. Note that this option is usually used with review.assembly after the number and ordering of chromosome-length scaffolds has been verified."

([ "$stage" != "" ] && ([ "$early" == "true" ] || [ "$fast" == "true" ])) && echo >&2 ":( Flags -e|--early--exit and -f|--fast-start are incompatible with setting a specific pipeline stage. Exiting!" && exit 1

############### HANDLE EXTERNAL DEPENDENCIES ###############

##	GNU Parallel Dependency
parallel="false"
if hash parallel 2>/dev/null; then
        ver=`parallel --version | awk 'NR==1{print \$3}'`
        [ $ver -ge 20150322 ] && parallel="true"
fi
[ $parallel == "false" ] && echo ":| WARNING: GNU Parallel version 20150322 or later not installed. We highly recommend to install it to increase performance. Starting pipeline without parallelization!" >&2

############### HANDLE ARGUMENTS ###############

[ -z $1 ] || [ -z $2 ] && echo >&2 "Not sure how to parse your input: files not listed or not found at expected locations. Exiting!" && echo >&2 "$USAGE_short" && exit 1

[ ! -s $1 ] || [ ! -s $2 ] && echo >&2 "Not sure how to parse your input: files not listed or not found at expected locations. Exiting!" && echo >&2 "$USAGE_short" && exit 1

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

if [ ! -f ${genomeid}.assembly ]; then
	awk -f ${pipeline}/utils/generate-assembly-file-from-fasta.awk ${orig_fasta} > ${genomeid}.assembly
fi

orig_assembly=${genomeid}.assembly

if grep -q ':::fragment_\|:::debris\|:::overhang_\|:::gap' "$orig_assembly"; then
	echo >&2 "The original fasta contains sequence labels that may interfere with labels assigned by the assembly pipeline. Please consider renaming your original sequences. Exiting!"
	exit 1
fi

awk -f ${pipeline}/utils/extract-cprops-from-assembly.awk ${orig_assembly}
orig_cprops=${genomeid}.cprops

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

############### ITERATIVE SCAFFOLDING/MISJOIN CORRECTION ###############

if [ "$stage" != "polish" ] && [ "$stage" != "split" ] && [ "$stage" != "seal" ] && [ "$stage" != "finalize" ]; then
	
	ROUND=0

	if ! [ "$fast" == "true" ]; then
        if [ -f ${genomeid}.*.cprops ] || [ -f ${genomeid}.mnd.*.txt ] ; then
            echo >&2 ":( Please remove or rename files ${genomeid}.*.cprops ${genomeid}.mnd.*.txt. Exiting!" && exit
        else
            ln -sf ${orig_cprops} ${genomeid}.${ROUND}.cprops
            ln -sf ${orig_mnd} ${genomeid}.mnd.${ROUND}.txt
        fi
	else
        ln -sf ${orig_mnd} ${genomeid}.mnd.${ROUND}.txt
		[ ! -f ${genomeid}.0.cprops ] || [ ! -f ${genomeid}.0.asm ] || [ ! -f ${genomeid}.0.hic ] || [ ! -f ${genomeid}.mnd.0.txt ] || [ ! -f ${genomeid}.0_asm.scaffold_track.txt ] || [ ! -f ${genomeid}.0_asm.superscaf_track.txt ] && echo >&2 ":( No early exit files are found. Please rerun the pipeline to include the round 0 assembly. Exiting!" && exit 1
	fi

	current_cprops=${genomeid}.${ROUND}.cprops
	current_mnd=${genomeid}.mnd.${ROUND}.txt

	echo "###############" >&1
	echo "Starting iterating scaffolding with editing:" >&1
	
	# run liger
	while true; do
	{
		if !([ "$fast" == "true" ] && [ ${ROUND} -eq 0 ] ); then
        # scaffold
            echo "...starting round ${ROUND} of scaffolding:" >&1
            bash ${pipeline}/scaffold/run-liger-scaffolder.sh -p ${parallel} -s ${input_size} -q ${mapq} ${current_cprops} ${current_mnd}
            
        # build a hic map of the resulting assembly
            echo "...visualizing round ${ROUND} results:" >&1
            bash ${pipeline}/visualize/run-asm-visualizer.sh -p ${parallel} -q ${mapq} -i -c ${current_cprops} ${genomeid}.${ROUND}.asm ${current_mnd}
			rm ${current_mnd}
		fi

		# break out of the scaffold-mismatch detection loop if the max number of steps is reached
			[ ${ROUND} -eq ${MAX_ROUNDS} ] && break

		# annotate near-diagonal mismatches in the map
			echo "...detecting misjoins in round ${ROUND} assembly:" >&1
			bash ${pipeline}/edit/run-mismatch-detector.sh -p ${parallel} -c ${editor_saturation_centile} -w ${editor_coarse_resolution} -d ${editor_coarse_region} -k ${editor_coarse_stringency} -n ${editor_fine_resolution} -b ${editor_coarse_norm} ${genomeid}.${ROUND}.hic
		# annotate repeats by coverage analysis
			bash ${pipeline}/edit/run-coverage-analyzer.sh -w ${editor_coarse_resolution} -t ${editor_repeat_coverage} ${genomeid}.${ROUND}.hic
		# store intermediate mismatch stuff	- not necessary
			mv depletion_score_wide.wig depletion_score_wide.at.step.${ROUND}.wig
			mv depletion_score_narrow.wig depletion_score_narrow.at.step.${ROUND}.wig
			mv mismatch_wide.bed mismatch_wide.at.step.${ROUND}.bed
			mv mismatch_narrow.bed mismatch_narrow.at.step.${ROUND}.bed
		# store intermediate repeat stuff - not necessary
			mv coverage_wide.wig coverage_wide.at.step.${ROUND}.wig
			mv repeats_wide.bed repeats_wide.at.step.${ROUND}.bed

			annotated=$(awk '{c+=$3-$2}END{print c}' mismatch_narrow.at.step.${ROUND}.bed)
			printf "...total %% genome annotated during mismatch analysis is: %.2f%%.\n" "$((10**2 * annotated*scale*100/totlength))e-2"
			if [ $((annotated*scale*100/totlength)) -gt 10 ]; then
				echo >&2 ":| WARNING: more than 10% of your genome assembly is flagged as problematic during mismatch analysis. Maybe your data is too sparse to do misjoin correction at default resolution. Please examine the relevant tracks and consider adjusting the --editor-coarse-resolution flag or other relevant editor flags!"
			fi

			annotated=$(awk '{c+=$3-$2}END{print c}' repeats_wide.at.step.${ROUND}.bed)
			printf "...total %% genome annotated during coverage analysis is %.2f%%.\n" "$((10**2 * annotated*scale*100/totlength))e-2"
			if [ $((annotated*scale*100/totlength)) -gt 10 ]; then
				echo >&2 ":| WARNING: more than 10% of your genome assembly is flagged as problematic during coverage analysis. Perhpaps your assembly contains a large number of collapsed repeats, alt heterozygosity or contaminan sequences. Please examine the relevant tracks and consider adjusting the --editor-repeat-coverage flag or other relevant editor flags!"
			fi

			resolved=$(awk 'NR==2{print $3}' ${genomeid}.${ROUND}_asm.superscaf_track.txt)  # scaled coordinates    
			awk -v end_interval=resolved -v bin=${editor_coarse_resolution} -f ${pipeline}/supp/plot_coverage.awk coverage_wide.at.step.${ROUND}.wig > coverage_wide.at.step.${ROUND}.dist.txt

		# consolidate bed annotations
			cat mismatch_narrow.at.step.${ROUND}.bed repeats_wide.at.step.${ROUND}.bed | sort -k 2,2n | awk 'BEGIN{FS="\t"; OFS="\t"}NR==1{start=$2; end=$3; next}$2<=end{if($3>end){end=$3}; next}{print "assembly", start, end; start=$2; end=$3}END{print "assembly", start, end}' > suspect.at.step.${ROUND}.bed

		# convert bed track into 2D annotations
			resolved=$(awk 'NR==2{print $3}' ${genomeid}.${ROUND}_asm.superscaf_track.txt)	# scaled coordinates	
			awk -v bin_size=${editor_fine_resolution} -f ${pipeline}/edit/overlay-edits.awk ${genomeid}.${ROUND}_asm.scaffold_track.txt suspect.at.step.${ROUND}.bed | awk -v r=${resolved} 'NR==1||$3<=r' > suspect_2D.at.step.${ROUND}.txt

		# separate intra and inter-input scaffold mismatches
			awk 'NR==1||$8=="debris"' suspect_2D.at.step.${ROUND}.txt > edits.for.step.$((ROUND+1)).txt
		# optional
			awk 'NR==1||$8=="mismatch"' suspect_2D.at.step.${ROUND}.txt > mismatches.at.step.$ROUND.txt

		# early exit on round zero if requested
			if [ "$early" == "true" ]; then
				echo >&1 ":) Early exit requested with --early-exit flag. Exiting after the first scaffolding round!"
				exit 0
			fi

		# check if there are any more edits to be done
			test=`wc -l < edits.for.step.$((ROUND+1)).txt`
			[ $test -eq 1 ] && echo >&1 ":) No more input edits to be done. Moving to polishing!" && rm edits.for.step.$((ROUND+1)).txt && break

		# move on to the next step
			ROUND=$((ROUND+1))
			[ -f ${genomeid}".edits.txt" ] && cp ${genomeid}".edits.txt" "archive."${genomeid}".edits.at.step."$((ROUND-1))".txt" # not necessary
		
		# reconstruct current edits
			awk 'BEGIN{OFS="\t"; print "chr1", "x1", "x2", "chr2", "y1", "y2", "color", "id", "X1", "X2", "Y1", "Y2"}$1~/:::debris/{print $1, 0, $3, $1, 0, $3, "0,0,0", "debris", 0, $3, 0, $3}' ${current_cprops} | awk -f ${pipeline}/lift/lift-input-annotations-to-asm-annotations.awk ${current_cprops} <(awk '{print $2}' ${current_cprops}) - | awk -f ${pipeline}/lift/lift-asm-annotations-to-input-annotations.awk ${orig_cprops} <(awk '{print $2}' ${orig_cprops}) - > h.old.edits.txt

		# add new edits
			bash ${pipeline}/lift/lift-edit-asm-annotations-to-original-input-annotations.sh ${orig_cprops} ${current_cprops} ${genomeid}.$((ROUND-1)).asm edits.for.step.${ROUND}.txt > h.new.edits.txt
			awk 'NR==1' "h.new.edits.txt" > temp		
			{ awk 'NR>1' h.old.edits.txt ; awk 'NR>1' "h.new.edits.txt" ; } | sort -k 1,1 -k 2,2n >> temp
			mv temp ${genomeid}".edits.txt"
			rm h.old.edits.txt h.new.edits.txt
		
		# apply edits		
			bash ${pipeline}/edit/apply-edits-prep-for-next-round.sh -p ${parallel} -r ${ROUND} ${genomeid}".edits.txt" ${orig_cprops} ${orig_mnd}
			current_cprops=${genomeid}.${ROUND}.cprops
			current_mnd=${genomeid}.mnd.${ROUND}.txt
	}
	done

	ln -sf ${genomeid}.${ROUND}.cprops ${genomeid}.resolved.cprops
	ln -sf ${genomeid}.${ROUND}.asm ${genomeid}.resolved.asm
	ln -sf ${genomeid}.${ROUND}_asm.scaffold_track.txt ${genomeid}.resolved_asm.scaffold_track.txt
	ln -sf ${genomeid}.${ROUND}_asm.superscaf_track.txt ${genomeid}.resolved_asm.superscaf_track.txt
	ln -sf ${genomeid}.${ROUND}.hic ${genomeid}.resolved.hic
	#ln -sf ${genomeid}.mnd.${ROUND}.txt ${genomeid}.mnd.resolved.txt

fi


############### POLISHING ###############

if [ "$stage" != "split" ] && [ "$stage" != "seal" ] && [ "$stage" != "merge" ] && [ "$stage" != "finalize" ]; then

	[ ! -f ${genomeid}.resolved.cprops ] || [ ! -f ${genomeid}.resolved.asm ] || [ ! -f ${genomeid}.resolved.hic ] && echo >&2 ":( No resolved files are found. Please rerun the pipeline to include the scaffold segment. Exiting!" && exit 1

	echo "###############" >&1
	echo "Starting polish:" >&1
	
	bash ${pipeline}/polish/run-asm-polisher.sh -p ${parallel} -q ${mapq} -j ${genomeid}.resolved.hic -a ${genomeid}.resolved_asm.scaffold_track.txt -b ${genomeid}.resolved_asm.superscaf_track.txt -s ${polisher_input_size} -c ${polisher_saturation_centile} -w ${polisher_coarse_resolution} -d ${polisher_coarse_region} -k ${polisher_coarse_stringency} -n ${polisher_fine_resolution} ${genomeid}.cprops ${orig_mnd} ${genomeid}.resolved.cprops ${genomeid}.resolved.asm
	
	mv ${genomeid}.resolved.polish.cprops ${genomeid}.polished.cprops
	mv ${genomeid}.resolved.polish.asm ${genomeid}.polished.asm

	mv ${genomeid}.resolved.polish.edits_2D.txt ${genomeid}.polished.edits_2D.txt 
	mv ${genomeid}.resolved.polish.mismatches_2D.txt ${genomeid}.polished.mismatches_2D.txt
	mv ${genomeid}.resolved.polish.suspect_2D.txt ${genomeid}.polished.suspect_2D.txt
	mv ${genomeid}.resolved.polish.mismatch_narrow.bed ${genomeid}.polished.mismatch_narrow.bed
	mv ${genomeid}.resolved.polish.depletion_score_narrow.wig ${genomeid}.polished.depletion_score_narrow.wig
	mv ${genomeid}.resolved.polish.mismatch_wide.bed ${genomeid}.polished.mismatch_wide.bed
	mv ${genomeid}.resolved.polish.depletion_score_wide.wig ${genomeid}.polished.depletion_score_wide.wig
	
	mv ${genomeid}.resolved.polish.hic ${genomeid}.polished.hic
	mv ${genomeid}.resolved.polish_asm.superscaf_track.txt ${genomeid}.polished_asm.superscaf_track.txt
	mv ${genomeid}.resolved.polish_asm.scaffold_track.txt ${genomeid}.polished_asm.scaffold_track.txt

fi

############### SPLITTING ###############

if [ "$stage" != "seal" ] && [ "$stage" != "merge" ] && [ "$stage" != "finalize" ]; then
	
	[ ! -s ${genomeid}.polished.cprops ] || [ ! -s ${genomeid}.polished.asm ] && echo >&2 ":( No resolved files are found. Please rerun the pipeline to include the scaffold/scaffold+polish segment. Exiting!" && exit 1
	
#	[ $chrom_num -ne 1 ] && bash ${pipeline}/split/run-asm-splitter.sh -c ${chrom_num} -r ${diploid} ${genomeid}.polished.cprops ${genomeid}.polished.asm ${genomeid}.mnd.polished.txt || cp ${genomeid}.polished.cprops ${genomeid}.polished.split.asm

	echo "###############" >&1
	echo "Starting split:" >&1
	bash ${pipeline}/split/run-asm-splitter.sh -p ${parallel} -q ${mapq} -j ${genomeid}.polished.hic -a ${genomeid}.polished_asm.scaffold_track.txt -b ${genomeid}.polished_asm.superscaf_track.txt -s ${splitter_input_size} -c ${splitter_saturation_centile} -w ${splitter_coarse_resolution} -d ${splitter_coarse_region} -k ${splitter_coarse_stringency} -n ${splitter_fine_resolution} ${genomeid}.cprops ${orig_mnd} ${genomeid}.polished.cprops ${genomeid}.polished.asm

	mv ${genomeid}.polished.split.cprops ${genomeid}.split.cprops
	mv ${genomeid}.polished.split.asm ${genomeid}.split.asm
	mv ${genomeid}.polished.split.hic ${genomeid}.split.hic
	mv ${genomeid}.polished.split_asm.superscaf_track.txt ${genomeid}.split_asm.superscaf_track.txt
	mv ${genomeid}.polished.split_asm.scaffold_track.txt ${genomeid}.split_asm.scaffold_track.txt
	
fi

############### SEALING ###############

if [ "$stage" != "merge" ] && [ "$stage" != "finalize" ]; then

	if [ ! -z ${review_assembly} ]; then
		ln -sf ${review_assembly} ${genomeid}.split.assembly
		awk -f ${pipeline}/utils/convert-assembly-to-cprops-and-asm.awk ${genomeid}.split.assembly
		unset review_assembly
	fi

	[ ! -s ${genomeid}.split.cprops ] || [ ! -s ${genomeid}.split.asm ] && echo >&2 ":( No split files are found. Please rerun the pipeline to include the split segment. Exiting!" && exit 1
	
	echo "###############" >&1
	echo "Starting sealing:" >&1
	
	# start slowly converting to .assembly as main input. If necessary split inside a particular block and cleanup after
	cat <(awk '{$0=">"$0}1' ${genomeid}.split.cprops) ${genomeid}.split.asm > ${genomeid}.split.assembly
	
	bash ${pipeline}/seal/run-assembly-sealer.sh -i ${input_size} ${genomeid}.split.assembly
	
	# prep for finalizing
	
	mv ${genomeid}.split.sealed.assembly ${genomeid}.rawchrom.assembly
	
fi

############### FINALIZING ###############

# finalize fasta

echo "###############" >&1
echo "Finalizing output:" >&1

if [ ! -z ${review_assembly} ] ; then
	cmp -s ${genomeid}.rawchrom.assembly $review_assembly || (rm -f ${genomeid}.rawchrom.assembly && ln -sf ${review_assembly} ${genomeid}.rawchrom.assembly)
fi

# sort assembly
awk -v n=${sort_chromosomes} -f $pipeline/utils/sort-assembly-by-size.awk ${genomeid}.rawchrom.assembly > ${genomeid}.rawchrom.assembly.tmp && mv ${genomeid}.rawchrom.assembly.tmp ${genomeid}.rawchrom.assembly
awk -f ${pipeline}/utils/convert-assembly-to-cprops-and-asm.awk ${genomeid}.rawchrom.assembly

# build final map if different from rawchrom
bash ${pipeline}/edit/edit-mnd-according-to-new-cprops.sh -q ${genomeid}.rawchrom.cprops ${orig_mnd} > ${genomeid}.rawchrom.mnd.txt
bash ${pipeline}/visualize/run-assembly-visualizer.sh -p ${parallel} -q ${mapq} -i -c ${genomeid}.rawchrom.assembly ${genomeid}.rawchrom.mnd.txt
rm ${genomeid}.rawchrom.mnd.txt

# build final fasta
# rely on first check to silence the usual annoying warnings 
awk -f ${pipeline}/utils/wrap-fasta-sequence.awk ${orig_fasta} | awk -v label1=":::fragment_" -v label2=":::debris" -f ${pipeline}/edit/edit-fasta-according-to-new-cprops.awk ${genomeid}.rawchrom.cprops - > ${genomeid}.rawchrom.fasta

options="-s ${input_size} -l ${genomeid} -g ${gap_size}"
[ ! -z ${chromosome_count} ] && options=$options" -c ${chromosome_count}"
[[ ! -z ${organism} ]] && options=$(printf "%s -o %q" "$options" "$organism")
[[ ! -z ${isolate} ]] && options=$(printf "%s -i %q" "$options" "$isolate")

eval "bash ${pipeline}/finalize/finalize-output.sh ${options} ${genomeid}.rawchrom.cprops ${genomeid}.rawchrom.asm ${genomeid}.rawchrom.fasta final"
rm ${genomeid}.rawchrom.fasta ${genomeid}*.cprops ${genomeid}*.asm

# if requested build _HiC.hic
if [ ! -z ${chromosome_count} ]; then
	bash ${pipeline}/visualize/build-sandboxed-hic.sh -c ${chromosome_count} ${genomeid}_HiC.assembly ${genomeid}.mnd.txt
fi
