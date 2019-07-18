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
path_to_input_fasta			Specify file path to draft assembly fasta file.
path_to_input_mnd			Specify path to deduplicated list of alignments of Hi-C reads to the draft assembly fasta as produced by the Juicer pipeline: the merged_nodups file (mnd).

OPTIONS:
-c|--chr chr
		Phase only one molecule (default phases all molecules listed in vcf file).
-q|--mapq mapq
		Consider only Hi-C reads that align with minimal mapping quality of mapq (default is 1).
-s|--stringency stringency
		Specify stringency parameter for the phaser (default is 3, i.e. 3-fold enrichment in Hi-C signal to one molecule as compared to the other is necessary to phase).
-r|--relaxation relaxation
		Specify relaxation parameter for the phaser (default is 1, i.e. consider enrichment above the noice level of 1 read)
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
relaxation=1
verbose=""

############### HANDLE OPTIONS ###############

while :; do
	case $1 in
		-h|--help)
			echo "$USAGE" >&1
			exit 0
        ;;
        -c|--chr) OPTARG=$2
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
        -r|--relaxation) OPTARG=$2
			re='^[0-9]+$'
			if [[ $OPTARG =~ $re ]]; then
				echo "... -r|--relaxation flag was triggered, phasing requiring enrichment against $OPTARG read(s) background noise." >&1
				relaxation=$OPTARG
			else
				echo ":( Wrong syntax for relaxation parameter. Exiting!" >&2
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

#1) parse original vcf file [ potentially unnecessary safeguard against unsorted vcf ]
echo ":) Parsing vcf file..." >&1
[ -z ${psf} ] && ( cat ${vcf} | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' | awk -v chr=${chr} -v output_prefix="in" -f ${pipeline}/phase/vcf-to-psf-and-assembly.awk && mv "in.assembly" `basename ${vcf} .vcf`".in.assembly" ) && psf="in.psf"
echo ":) Done parsing vcf file." >&1

#2) extract SNP overlapping edges from mnd file
echo ":) Extracting Hi-C contacts overlapping SNPs...">&1
[ -z ${edge_mnd} ] && { bash ${pipeline}/phase/extract-SNP-edges-from-mnd-file.sh -q ${mapq} -o "snp.mnd.txt" ${psf} ${mnd} 2>&3 | sed 's/^/.../';  } 3>&1 1>&2 | sed 's/^/.../' && edge_mnd="snp.mnd.txt"
echo ":) Done extracting Hi-C contacts overlapping SNPs.">&1

#3) visualize input
echo ":) Visualizing input phased blocks..."
[ -f `basename ${vcf} .vcf`".in.assembly" ] && { bash ${pipeline}/visualize/run-assembly-visualizer.sh -c -p ${parallel} `basename ${vcf} .vcf`".in.assembly" ${edge_mnd} 2>&3 | sed 's/^/.../';  } 3>&1 1>&2 | sed 's/^/.../'
echo ":) Done visualizing input phased blocks."

#4) phase
echo ":) Phasing..."
if [ ! -z "$chr" ]; then
	{ awk -v stringency=${stringency} -v relaxation=${relaxation} -v outfile="out.psf" -v verbose=${verbose} -f ${pipeline}/phase/phase-intrachromosomal.awk ${psf} ${edge_mnd} 2>&3 | sed 's/^/.../';  } 3>&1 1>&2 | sed 's/^/.../'
else
	export SHELL=$(type -p bash)
	export psf=${psf}
	export edge_mnd=${edge_mnd}
	export stringency=${stringency}
	export relaxation=${relaxation}
	export verbose=${verbose}
	export pipeline=${pipeline}
	doit () { 
		cmd="awk -v chr=$1 '\$1==\">\"chr{print; id[\$NF]=1; id[-\$NF]=1}\$1~/^>/{next}(\$1 in id){print}' ${psf} > h.$1.psf && awk -v stringency=${stringency} -v relaxation=${relaxation} -v outfile=out.$1.psf -v verbose=${verbose} -f ${pipeline}/phase/phase-intrachromosomal.awk h.$1.psf ${edge_mnd} && rm h.$1.psf"
		eval $cmd
	}
	export -f doit
	
	if [ $parallel == "true" ]; then
		awk '$0~/^>/{print substr($1,2)}' ${psf} | sort -u | parallel --will-cite doit
		awk '$0~/^>/{print substr($1,2)}' ${psf} | sort -u | parallel --will-cite -k "awk '\$0~/^>/' out.{}.psf" > out.psf
		awk '$0~/^>/{print substr($1,2)}' ${psf} | sort -u | parallel --will-cite -k "awk '\$0!~/^>/' out.{}.psf" >> out.psf
		awk '$0~/^>/{print substr($1,2)}' ${psf} | sort -u | parallel --will-cite rm out.{}.psf
	else
		rm -f "out.psf.p1" "out.psf.p2"
		while read -r var; do
			doit ${var}
			awk '$0~/^>/' "out."${var}".psf" >> out.psf.p1
			awk '$0!~/^>/' "out."${var}".psf" >> out.psf.p2
			rm out.${var}.psf
		done < <(awk '$0~/^>/{print substr($1,2)}' ${psf} | sort -u)
		cat "out.psf.p1" "out.psf.p2" > "out.psf"
		rm "out.psf.p1" "out.psf.p2"
	fi
fi
echo ":) Done phasing."

#5) visualize results and dump vcf
echo ":) Visualizing output phased blocks..."
awk -f ${pipeline}/phase/psf-to-assembly.awk "out.psf" > `basename ${vcf} .vcf`".out.assembly" 
{ bash ${pipeline}/visualize/run-assembly-visualizer.sh -c -p ${parallel} `basename ${vcf} .vcf`".out.assembly" ${edge_mnd} 2>&3 | sed 's/^/.../';  } 3>&1 1>&2 | sed 's/^/.../'
awk -f ${pipeline}/phase/psf-to-vcf.awk "out.psf" > `basename ${vcf} .vcf`".out.vcf"
echo ":) Done visualizing output phased blocks."

#6) cleanup
#rm ${psf} ${edge_mnd} "out.psf"

