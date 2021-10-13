#!/bin/bash
## bash lift-input-mnd-to-HiC-mnd.sh <assembly> <mnd>
## Handle options

USAGE="
**********************

A parallelized merger of various script to edit and lift merged_nodups file (long format) according to an assembly file.

	./lift-mnd-file.sh [-q mapq] [-j jobcount] [-a] [-l] [-h] <path_to_assembly_file> <path_to_mnd_txt>

Input: assembly file and mnd file.

Options:
    -j
                specified number of jobs for parallelization. Acceptible formats: 2, \"80%\" (see GNU parallel for more details).
	-q  <mapq>:
                ignore reads with mapping quality less than mapq. Default: 1. Note that by default this means some reads will be lost (mapq=0).
    -a
                use \"assembly\" output format for output chromosome name and coordinate (i.e. genome-wide for \"assembly\" chromosome).
    -l
                keep to long format for output mnd. [Default: short format is used.]
    -c
                reports edges combinatorially in cases when single reads overlap multiple SNPs. 
    -o  <outfilename>
                output filename
    -h
                shows this help

************************
"

jobcount="80%"
sandbox=1
longformat=0
mapq=1
scale=1
outfilename="lifted.mnd.txt"

while getopts "j:q:alh" opt; do
case $opt in
    j)  #TODO: check if compatible with GNU parallel syntax
        echo "... -j flag was triggered, will try to pass -j $OPTARG to GNU parallel to set the number of simultaneous jobs." >&1
        jobs=$OPTARG
    ;;
    q)  re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]]; then
            echo "... -q flag was triggered, will ingore all reads in the mnd file with mapping quality lower than $OPTARG on either of the reads." >&1
            mapq=$OPTARG
        else
            echo ":( Wrong syntax for mapping quality. Exiting!" >&2
            exit 1
        fi
    ;;
    a)
        echo "... -a flag was triggered, will lift onto an \"assembly\" chromosome."
        sandbox=0
    ;;
    l)
        echo "... -l flag was triggered, will keep optional fields in the mnd file."
        longformat=1
    ;;
    o)
        echo "... -o flag was triggered, will use $OPTARG for output filename."
        outfilename=$OPTARG
    ;;
    h)
        echo "$USAGE" >&1
        exit 0
    ;;
    *) echo "$USAGE" >&2
        exit 1
    ;;
esac
done

shift $(( OPTIND-1 ))

## HANDLE ARGUMENTS

([ -z ${1} ] || [ -z ${2} ]) && (":( Some input is missing. Exiting!" >&2 && exit 1)

assembly=$1
mnd=$2

## DEPENDENCIES

pipeline=`cd "$( dirname $0)" && cd .. && pwd`
lift_mnd_script=${pipeline}/lift/lift-input-mnd-to-HiC-mnd.awk
##	GNU Parallel Dependency
parallel="false"
if hash parallel 2>/dev/null; then
        ver=`parallel --version | awk 'NR==1{print \$3}'`
        [ $ver -ge 20150322 ] && parallel="true"
fi

[ $parallel == "false" ] && echo >&2 ":| WARNING: GNU Parallel version 20150322 or later not installed. We highly recommend to install it to increase performance. Starting pipeline without parallelization!"

## Handle default zoom if not specified

if [ $sandbox -eq 0 ]; then
	# calculate necessary zoom
	totlength=`awk '$0~/^>/{len[$2]=$3;len[-$2]=$3;next}{for(i=1;i<=NF;i++){total+=len[$i]}}END{print total}' ${assembly}`
	scale=$(( 1 + $totlength / 2100000000 ))
    if [ $scale -ne 1 ]; then
        echo ":| Warning: scaling assembly to fit into Java int. This is expected for genomes > 2.1Gb."
    fi
fi

optstring="-v scale=${scale} -v sandbox=${sandbox} -v longformat=${longformat} -v mapq=${mapq}"

[ $sandbox -eq 1 ] && optstring=${optstring}" -v chromlabel=\"null\""

[ $parallel == "true" ] && parallel -a ${mnd} --pipepart --will-cite --jobs ${jobcount} --block 1G "awk ${optstring} -f ${lift_mnd_script} ${assembly} - " > ${outfilename} || awk ${optstring} -f ${lift_mnd_script} ${assembly} ${mnd} > ${outfilename}

[ $sandbox -eq 1 ] && sort -k2,2n -k6,6n --parallel=48 -S32G ${outfilename} | awk '{$2="HiC_scaffold_"$2; $6="HiC_scaffold_"$6}1' > ${outfilename}".tmp" && mv ${outfilename}".tmp" ${outfilename}