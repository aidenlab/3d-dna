#!/bin/bash

#### Description: Script to build aggregate chromosome contact maps.
#### Usage: bash ./build-aca-hic.sh [options] <path_to_chrom_sizes_file> <path_to_centromere_bed_or_bedpe> <path_to_mnd_file>.
#### Input: chrom.sizes file, file with centromere positions in bed or bedpe format and merged_nodups file. For the mnd file long format and short with score format are acceptable.
#### Output: hic file.
#### Parameters: mapq threshold (default is 1 for mapq>=1).
#### Dependencies: Java; GNU Parallel; 3D-DNA scripts.
#### Written by: Olga Dudchenko, version date 07/29/2020

USAGE="
*****************************************************
Building aggregate chromosome contact maps: 5/24/20.

USAGE: ./build-aca-hic.sh [options] <path_to_chrom_sizes_file> <path_to_centromere_bed_or_bedpe> <path_to_mnd_file>

DESCRIPTION:
This is a script to generate ACA type hic contact files.

ARGUMENTS:
path_to_chrom_sizes_file
						Specify path to chrom.sizes file describing which sequences among those appearing in the mnd are chromosomes and what is their length.
path_to_centromere_bed_or_bedpe
						Specify path to file annotating the position of chromosomes in the centromeres. (It is ok if these are listed only for a subset of chromosomes. Those without will be ignored.).

OPTIONS:
-q mapq
						Build map for a specific mapq threshold (default is 1). 
-h
						Shows this help.
*****************************************************
"
## 3D-DNA
pipeline=`cd "$( dirname $0)" && cd .. && pwd`

## Defaults
mapq=1

## HANDLE OPTIONS
while getopts "q:h" opt; do
case $opt in
	h) echo "$USAGE"
		exit 0
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
	*)  echo ":( Illegal options. Exiting."
		echo "$USAGE"
		exit 1
	;;
esac
done

shift $(( OPTIND-1 ))

## HANDLE ARGUMENTS
chrom_sizes_file=$1
centromere_bed_file=$2
merged_nodups_file=$3

export centromere_bed_file
export chrom_sizes_file
export mapq

doit() {
	awk -v mapq=${mapq} -f <(cat - <<-'EOD'
# read in chrom.sizes
FILENAME==ARGV[1]{
len[$1]=$2
next
}
# read centromeric bed/bedpe
FILENAME==ARGV[2]&&$0~/^#/{
next
}
FILENAME==ARGV[2]{
#filter for metacentric
if(($2+$3)>2/3*len[$1]&&($2+$3)<4/3*len[$1])
	cen[$1]=($2+$3)/2
next
}
# merged_nodups
!($2 in len) || !($6 in len) || !($2 in cen) || !($6 in cen) {next}
#mapq filtering if this is long format
NF>9&&($9<mapq||$12<mapq){next}
{
# record for rescale
A=cen[$2]
B=len[$2]-cen[$2]
C=cen[$6]
D=len[$6]-cen[$6]

#simplify

$1=0
$5=0

# short format vs long format
if(NF==9){score=$9}else{score=1}
$9=""
if($10){$10=""}
if($11){$11=""}
if($12){$12=""}
if($13){$13=""}
if($14){$14=""}
if($15){$15=""}
if($16){$16=""}

# handle intrachromosomal entries
if($2==$6)
{

$2="assembly"
$6="assembly"

if($3<=A)
$3=int($3/A*500000)
else
$3=500000+int(500000*($3-A)/B)

sym3=1000000-$3

if($7<=C)
$7=int($7/C*500000)
else
$7=500000+int(500000*($7-C)/D)

sym7=1000000-$7

print $0" "(score*2*(length(cen)-1))

$3+=1000000
$7+=1000000

print $0" "(score*2*(length(cen)-1))

$3=sym3
$7=sym7

print $0" "(score*2*(length(cen)-1))

$3+=1000000
$7+=1000000

print $0" "(score*2*(length(cen)-1))

}
else # handle interchromosomal entries
{

$2="assembly"
$6="assembly"
$4=0
$8=1

# control for more inter than intra

# move everything into the same quadrangle
if($3>A)
{
$3=A+B-$3
tmp=A
A=B
B=tmp
}
if($7>C)
{
$7=C+D-$7
tmp=C
C=D
D=tmp
}

# print 8 entries to equilibrate coverage biases

opos1=$3
opos2=$7

$3=int(opos1/A*500000)
$7=1000000+int(opos2/C*500000)
print $0" "score
$3=int(opos2/C*500000)
$7=1000000+int(opos1/A*500000)
print $0" "score


$3=int(opos1/A*500000)
$7=2000000-int(opos2/C*500000)
print $0" "score
$3=int(opos2/C*500000)
$7=2000000-int(opos1/A*500000)
print $0" "score

$3=1000000-int(opos1/A*500000)
$7=1000000+int(opos2/C*500000)
print $0" "score
$3=1000000-int(opos2/C*500000)
$7=1000000+int(opos1/A*500000)
print $0" "score

$3=1000000-int(opos1/A*500000)
$7=2000000-int(opos2/C*500000)
print $0" "score
$3=1000000-int(opos2/C*500000)
$7=2000000-int(opos1/A*500000)
print $0" "score
}
}
EOD
) $1 $2 $3
}

export -f doit

parallel -a ${merged_nodups_file} --pipepart -j 80% --will-cite --block 1G --env doit "doit ${chrom_sizes_file} ${centromere_bed_file} -" > interchr.mnd.txt

bash ${pipeline}/visualize/juicebox_tools.sh pre -q ${mapq} interchr.mnd.txt `basename ${chrom_sizes_file} .chrom.sizes`".aca.hic" <(echo "assembly	2000000")
