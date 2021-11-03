#!/bin/bash
##########
#The MIT License (MIT)
#
# Copyright (c) 2019 Aiden Lab
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
# Written by O.D. (olga.dudchenko@.bcm.edu)
# Phaser module for 3D-DNA pipeline
#

version=190716
echo "*****************************************************" >&1
echo "cmd log: "$0" "$* >&1
echo "*****************************************************" >&1

# 
# echo "version: "${version}

USAGE="
*****************************************************

USAGE: ./run-hic-phaser.sh [options] <path_to_vcf_file> <path_to_mnd_file> 

DESCRIPTION:
This is a script to assemble draft assemblies (represented in input by draft fasta and deduplicated list of alignments of Hi-C reads to this fasta as produced by the Juicer pipeline) into chromosome-length scaffolds. The script will produce an output fasta file, a Hi-C map of the final assembly, and a few supplementary annotation files to help review the result in Juicebox.

ARGUMENTS:
path_to_input_vcf			Specify file path to vcf file with unphased or partially phased SNPs.
path_to_input_mnd			Specify path to deduplicated list of alignments of Hi-C reads to the draft assembly fasta as produced by the Juicer pipeline: the merged_nodups file (mnd).

OPTIONS:
-c|--chr chr
		Phase only one molecule (default phases all molecules listed in vcf file).
-q|--mapq mapq
		Consider only Hi-C reads that align with minimal mapping quality of mapq (default is 1).
-s|--stringency stringency
		Specify stringency parameter for the phaser (default is 3, i.e. 3-fold enrichment in Hi-C signal to one molecule as compared to the other is necessary to phase).
-b|--background background
		Specify background parameter for the phaser (default is 1, i.e. calculate enrichment on top of the noise level of 1 read)
-h|--help
		Shows this help.


TODO: need to pass filteres to parcing vcf file. Currently: all \"PASS\" variants are considered.
COULDDO: pass a name for multi-sample vcf files.

*****************************************************
"

# path to 3D-DNA pipeline
pipeline=`cd "$( dirname $0)" && cd .. && pwd`

# defaults:
chr=""
mapq=1
stringency=3
background=1
verbose=""
build_maps_in_native_coordinates=false

############### HANDLE OPTIONS ###############

while :; do
	case $1 in
		-h|--help)
			echo "$USAGE" >&1
			exit 0
        ;;
        -c|--chr) OPTARG=$2
			echo "... -c|--chr flag was triggered, ignoring all sequences in the vcf except for $OPTARG." >&1
			chr=$OPTARG
        	shift
        ;;
        -q|--mapq) OPTARG=$2
			re='^[0-9]+$'
			if [[ $OPTARG =~ $re ]]; then
				echo "... -q|--mapq flag was triggered, phasing using reads with at least $OPTARG mapping quality." >&1
				mapq=$OPTARG
			else
				echo ":( Wrong syntax for mapping quality parameter. Exiting!" >&2
				exit 1
			fi
        	shift
        ;;
        -s|--stringency) OPTARG=$2
			re='^[0-9]+$'
			if [[ $OPTARG =~ $re ]]; then
				echo "... -s|--stringency flag was triggered, phasing requiring at least $OPTARG enrichment of one haplotype vs the other." >&1
				stringency=$OPTARG
			else
				echo ":( Wrong syntax for stringency parameter. Exiting!" >&2
				exit 1
			fi
        	shift
        ;;
        -b|--background) OPTARG=$2
			re='^[0-9]+$'
			if [[ $OPTARG =~ $re ]]; then
				echo "... -b|--background flag was triggered, phasing requiring enrichment against $OPTARG read(s) background noise." >&1
				background=$OPTARG
			else
				echo ":( Wrong syntax for background parameter. Exiting!" >&2
				exit 1
			fi
        	shift
        ;;
### unprompted
		-p|--psf) OPTARG=$2
			psf=$OPTARG
			shift
		;;
		-e|--edges) OPTARG=$2
			edge_mnd=$OPTARG
			shift
		;;
		-v|--verbose) OPTARG=$2
			verbose=$OPTARG
			shift
		;;
		# --build-maps-in-native-coordinates)
		# 	build_maps_in_native_coordinates=true;
		# ;;
### utilitarian
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

[ -z $1 ] || [ -z $2 ] && echo >&2 ":( Not sure how to parse your input: files not listed or not found at expected locations. Exiting!" && echo >&2 "$USAGE" && exit 1

[ ! -s $1 ] || [ ! -s $2 ] && echo >&2 ":( Not sure how to parse your input: files not listed or not found at expected locations. Exiting!" && echo >&2 "$USAGE" && exit 1

if [ "$#" -ne 2 ]; then
    echo >&2 "Illegal number of arguments. Please double check your input. Exiting!" && echo >&2 "$USAGE" && exit 1
fi

vcf=$1
mnd=$2

############### MAIN #################

#1) parse original vcf file [ assumes the vcf file is sorted ]
echo ":) Parsing vcf file..." >&1
[ -z ${psf} ] && ( awk -v chr=${chr} -v output_prefix="in" -f ${pipeline}/phase/vcf-to-psf-and-assembly.awk ${vcf} && mv "in.assembly" `basename ${vcf} .vcf`".in.assembly" ) && psf="in.psf"
if grep -q ':' ${psf}; then
	echo >&2 ":( Some of the reference sequences appear to have a \":\" character in their names. This is going to interfere with the hic phaser internal annotation schema. Please rename the sequences in the vcf and the mnd files to proceed."
	echo ":( Something went wrong. Check stderr for more info. Exiting!" && exit 1
fi

echo ":) Done parsing vcf file." >&1

#2) extract SNP overlapping edges from mnd file
echo ":) Extracting Hi-C contacts overlapping SNPs...">&1
#[ -z ${edge_mnd} ] && { bash ${pipeline}/phase/extract-SNP-edges-from-mnd-file.sh -q ${mapq} -o "snp.mnd.txt" ${psf} ${mnd} 2>&3 | sed 's/^/.../';  } 3>&1 1>&2 | sed 's/^/.../' && edge_mnd="snp.mnd.txt"
if [ -z ${edge_mnd} ]; then
	{ bash ${pipeline}/phase/extract-SNP-reads-from-mnd-file.sh -q ${mapq} -d -n -o "dangling.native.mnd.txt" ${psf} ${mnd} | sed 's/^/.../'; } 2> >(while read line; do echo "...$line" >&2; done)
	[ ${PIPESTATUS[0]} -ne 0 ] && echo ":( Something went wrong. Check stderr for more info. Exiting!" && exit 1
#	awk '($2~/:/)&&($6~/:/)&&($2!=$6){$3=1;$7=1;print}' "dangling.native.mnd.txt" > "snp.mnd.txt"; 

		awk -v mapq=$mapq '($2~/:/)&&($6~/:/)&&($2!=$6)&&$9>=mapq&&$12>=mapq{$3=1;$7=1;print}' "dangling.native.mnd.txt" > "snp.mnd.txt"; 

	edge_mnd="snp.mnd.txt"
fi

echo ":) Done extracting Hi-C contacts overlapping SNPs.">&1

#3) visualize input
echo ":) Visualizing input phased blocks..."
# if [ -f `basename ${vcf} .vcf`".in.assembly" ]; then
# 	{ bash ${pipeline}/visualize/run-assembly-visualizer.sh -c -p ${parallel} `basename ${vcf} .vcf`".in.assembly" ${edge_mnd} | sed 's/^/.../'; } 2> >(while read line; do echo "...$line" >&2; done)
# 	[ ${PIPESTATUS[0]} -ne 0 ] && echo ":( Something went wrong. Check stderr for more info. Exiting!" && exit 1
# fi

echo ":) Done visualizing input phased blocks."

#4) phase
echo ":) Phasing..."
if [ ! -z "$chr" ] && [[ $chr != *"|"* ]]; then
	{ awk -v stringency=${stringency} -v background=${background} -v outfile="out.psf" -v verbose=${verbose} -f ${pipeline}/phase/phase-intrachromosomal.awk ${psf} ${edge_mnd} | sed 's/^/.../';  } 2> >(while read line; do echo "...$line" >&2; done)
	[ ${PIPESTATUS[0]} -ne 0 ] && echo ":( Something went wrong. Check stderr for more info. Exiting!" && exit 1
else
	if [ -z "$chr" ]; then
		chr=$(awk '$0~/^>/{if($1!=prev){str=str"|"substr($1,2); prev=$1;}}END{print substr(str,2)}' ${psf})
	fi
	
	export SHELL=$(type -p bash)
	export psf=${psf}
	export edge_mnd=${edge_mnd}
	export stringency=${stringency}
	export background=${background}
	export verbose=${verbose}
	export pipeline=${pipeline}
	doit () { 
		cmd="echo \"Phasing chr $1.\" && awk -v chr=$1 '\$1==\">\"chr{print; id[\$NF]=1; id[-\$NF]=1}\$1~/^>/{next}(\$1 in id){print}' ${psf} > h.$1.psf && awk -v stringency=${stringency} -v background=${background} -v outfile=out.$1.psf -v verbose=${verbose} -f ${pipeline}/phase/phase-intrachromosomal.awk h.$1.psf ${edge_mnd} && rm h.$1.psf"
		eval $cmd
	}
	export -f doit

	if [ $parallel == "true" ]; then
		echo $chr | tr "|" "\n" | parallel --will-cite doit
		echo $chr | tr "|" "\n" | parallel --will-cite -k "awk '\$0~/^>/' out.{}.psf" > out.psf
		echo $chr | tr "|" "\n" | parallel --will-cite -k "awk '\$0!~/^>/' out.{}.psf" >> out.psf
		echo $chr | tr "|" "\n" | parallel --will-cite rm out.{}.psf
	else
		rm -f "out.psf.p1" "out.psf.p2"
		while read -r var; do
			doit ${var}
			awk '$0~/^>/' "out."${var}".psf" >> out.psf.p1
			awk '$0!~/^>/' "out."${var}".psf" >> out.psf.p2
			rm out.${var}.psf
		done < <(echo $chr | tr "|" "\n")
		cat "out.psf.p1" "out.psf.p2" > "out.psf"
		rm "out.psf.p1" "out.psf.p2"
	fi
fi

#TODO: exit code tracking?
echo ":) Done phasing."

#5) visualize results and dump vcf
echo ":) Visualizing output phased blocks..."
awk -f ${pipeline}/phase/psf-to-assembly.awk "out.psf" > `basename ${vcf} .vcf`".out.assembly" 

if [ -f `basename ${vcf} .vcf`".out.assembly" ]; then
	{ bash ${pipeline}/visualize/run-assembly-visualizer.sh -c -p ${parallel} `basename ${vcf} .vcf`".out.assembly" ${edge_mnd} | sed 's/^/.../'; } 2> >(while read line; do echo "...$line" >&2; done)
	[ ${PIPESTATUS[0]} -ne 0 ] && echo ":( Something went wrong. Check stderr for more info. Exiting!" && exit 1
fi
echo ":) Done visualizing output phased blocks."

#6) if requested build contact maps for the largest component in 'native' coordinates. Optionally keep dangling ends.

bash ${pipeline}/phase/replace-variants-with-homologs-in-var-mnd.sh -k out.psf dangling.native.mnd.txt
awk -f ${pipeline}/phase/vcf-to-native-assembly.awk ${vcf} > native.assembly
bash ${pipeline}/lift/lift-input-mnd-to-HiC-mnd.sh native.assembly homolog.mnd.txt


#7) update the vcf file
awk -f ${pipeline}/phase/psf-to-vcf.awk "out.psf" > `basename ${vcf} .vcf`".out.vcf"

#8) cleanup
#rm ${psf} ${edge_mnd} "out.psf"

