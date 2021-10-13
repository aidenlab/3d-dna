#!/bin/bash
## Helper script to generate rainbow tracks
## Written by: OD olga.dudchenko@bcm.edu

USAGE="
*****************************
	./make-rainbow-tracks.sh -b <bin> [chr:start-end]
*****************************
"
bin=500000;
chr="assembly"
start=0;
end=2100000000;

## HANDLE OPTIONS

while getopts "hb:" opt; do
case $opt in
    h) echo "$USAGE" >&1
        exit 0
    ;;
	b) bin=$OPTARG
    ;;
    *) echo "$USAGE" >&2
        exit 1
    ;;
esac
done

shift $(( OPTIND-1 ))

## HANDLE ARGUMENTS

[ -z $1 ] || { chr=${1%%:*}; str=${1##*:}; start=${str%%-*}; end=${str##*-}; [ -z $chr ] || [ -z $start ] || [ -z $end ] && exit 1; }
size=$(($((end-start))/bin))

function sin_to_rgb {
    i=$1;
    PI=3.14159265
    IFS='   ' read -a arr <<<$(echo "var1=127*s($PI/$size*2*$i+0*$PI*2/3)+128; var2=127*s($PI/$size*2*$i+2*$PI*2/3)+128; var3=127*s($PI/$size*2*$i+1*$PI*2/3)+128; scale=0; var1/1; var2/1; var3/1" | bc -l)
    echo "${arr[0]},${arr[1]},${arr[2]}"
}

echo "track name=\"ItemRGBDemo\" description=\"Item RGB demonstration\" visibility=2 itemRgb=\"On\""
for i in $(seq 0 $size)
do
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $chr $((i*bin+start)) $(((i+1)*bin+start)) "-" 0 "+" $((i*bin+start)) $(((i+1)*bin+start)) `sin_to_rgb $i`
done



