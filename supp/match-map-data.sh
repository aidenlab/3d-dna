#!/bin/bash

USAGE="
*****************************
	./generate-linkage-map.sh -d <Juneja/...> -c <chrom_number> <path_to_orig_cprops> <path_to_final_tiled_cprops> <path_to_final_tiled_asm>
*****************************
"

## HANDLE OPTIONS

while getopts "hd:c:" opt; do
case $opt in
    h) echo "$USAGE" >&1
        exit 0
    ;;
    d) data=$OPTARG
    ;;
    c) re='^[0-9]+$'
    if [[ $OPTARG =~ $re ]] && [[ $OPTARG -gt 0 ]]; then
	chr=$OPTARG
    else
	echo ":( Wrong syntax for chromosome number. Exiting!" >&2 && exit 1
    fi
    ;;
    *) echo "$USAGE" >&2
        exit 1
    ;;
esac
done

shift $(( OPTIND-1 ))

[ -z $data ] && echo >&1 ":| Warning: Map data reference was not listed. Running against default data from Juneja et al, 2014. To request a different map type -d Juneja/..." && data=Juneja

pipeline=`cd "$( dirname $0)" && cd .. && pwd`


case "$data" in
	Juneja)
		datafile="$pipeline/data/Juneja_et_al_linkage_map.dat"
		chr=3
	;;
	Timoshevskiy)
		datafile="$pipeline/data/Timoshevskiy_et_al.dat"
		chr=3
	;;
	Hickner)
		datafile="$pipeline/data/Hickner_et_al_linkage_map.dat"
		chr=3
	;;
	Aresburger)
		datafile="$pipeline/data/Aresburger_et_al.dat"
		chr=3
	;;
	Unger)
		datafile="$pipeline/data/Unger_et_al.dat"
		chr=3
	;;
	Naumenko)
		datafile="$pipeline/data/Naumenko_et_al.dat"
		chr=3
	;;
	*) 
		datafile=$data
		[ -z $chr ] && echo "Unknown data or no chromosome number. Exiting!" >&2 && exit 1
	;;
esac

## HANDLE ARGUMENTS

[ -z $1 ] || [ -z $2 ] || [ -z $3 ] && echo >&2 ":( Some input seems missing." && echo >&2 "$USAGE" && exit 1

orig_cprops=$1
current_cprops=$2
current_asm=$3

id=`basename ${current_asm} .asm`

## TODO: check compatibility


awk 'BEGIN{FS="\t"; OFS="\t"; print "chr1", "x1", "x2", "chr2", "y1", "y2", "color", "id", "X1", "X2", "Y1", "Y2"}NR>1{print $4, $5-1, $6, $4, $5-1, $6, "0,0,0", "+"$1, $5-1, $6, $5-1, $6}' ${datafile} | awk -f ${pipeline}/lift/lift-input-annotations-to-asm-annotations.awk ${orig_cprops} <(awk '{print $2}' ${orig_cprops}) - | awk -f ${pipeline}/supp/lift-asm-annotations-to-input-annotations.awk ${current_cprops} <(awk '{print $2}' ${current_cprops}) - | awk -f ${pipeline}/lift/lift-input-annotations-to-asm-annotations.awk ${current_cprops} ${current_asm} - > ${data}_vs_${id}_2D_annotations.txt

awk -v chr=${chr} 'FILENAME==ARGV[1]{len[$2]=$3; next}FNR<=chr{gsub("-",""); c=0; for(i=1;i<=NF;i++){c+=len[$i]}; print FNR"\t"c}' ${current_cprops} ${current_asm} > ${id}.chrom.tiled.sizes

awk -v chr=${chr} 'BEGIN{FS="\t"; OFS="\t"}FILENAME==ARGV[1]{shift[$1+1]=shift[$1]+$2; next}FILENAME==ARGV[2]{pos[substr($8,2)]=$2+1; next}FNR==1{print $0, "Chr-length scaffold (current study)", "Start position (current study)"; next}{if (!pos[$1]){chrscaf="NA"} else {chrscaf=1; while(chrscaf<=chr){if(pos[$1]<shift[chrscaf+1]){break};chrscaf++}; if(chrscaf==chr+1){chrscaf="NA"}}; if(chrscaf=="NA"){print $0, "NA", "NA"}else{print $0, chrscaf, pos[$1]-shift[chrscaf], pos[$1]}}' ${id}.chrom.tiled.sizes $data"_vs_"${id}"_2D_annotations.txt" ${datafile} > ${data}_map_vs_${id}.txt
