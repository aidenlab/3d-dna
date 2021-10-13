#!/bin/bash
## This is a helper script to generate .dotplot.txt files from maf annotations 
## Usage: ./lift-maf-annotations-to-assembly-dotplot.sh [ -t number_of_scaffolds_to_include_for_target ] [ -q number_of_scaffolds_to_include_for_query ] <target.assembly> <query.assembly> <target_annotations> <query_annotations>
## Written by: Olga Dudchenko

USAGE=" ./lift-maf-annotations-to-assembly-dotplot.sh [ -t <first_t_scaffolds_in_target> ] [ -q <first_q_scaffolds_in_query> ] <target_assembly_file> <query_assembly_file> <maf_file> "

# HANDLE DEPENDENCIES

## 3D-DNA
pipeline=`cd "$( dirname $0)" && cd .. && pwd`

# HANDLE OPTIONS

bundle_size=15000

while getopts "t:q:b:h" opt; do
case $opt in
	t) target_number=$OPTARG
	;;
	q) query_number=$OPTARG
	;;
	h) echo "$USAGE"
		exit 0
	;;
	*) echo "$USAGE" >&2
		exit 1
	;;
esac
done
shift $(( OPTIND-1 ))

# HANDLE ARGUMENTS

target_assembly=$1
query_assembly=$2
maf=$3

awk -v target_file_name="target_annotations.txt" -v query_file_name="query_annotations.txt" -f ${pipeline}/supp/maf-to-annotations.awk ${maf}

target_annotations="target_annotations.txt"
query_annotations="query_annotations.txt"

## safeguards
([ -f ${target_assembly} ] && [ -f ${query_assembly} ] && [ -f ${target_annotations} ] && [ -f ${query_annotations} ] ) || (echo ":( Expected arguments not found. Exiting!" && exit 1)

[ -z ${target_number} ] && target_number=`awk '$0!~/^>/{counter++}END{print counter}' ${target_assembly}`
[ -z ${query_number} ] && query_number=`awk '$0!~/^>/{counter++}END{print counter}' ${query_assembly}`

# target
bash ${pipeline}/supp/edit-and-lift-input-annotations-to-assembly-annotations.sh -b ${bundle_size} ${target_assembly} ${target_annotations} | awk -f ${pipeline}/lift/lift-assembly-annotations-to-input-annotations.awk <(awk '$0~/^>/{len[$2]=$3;len[-$2]=$3;next}{counter++; for(i=1;i<=NF;i++){clen[counter]+=len[$i]}}END{for(i=1;i<=counter;i++){print ">"i, i, clen[i]};for(i=1;i<=counter;i++){print i}}' ${target_assembly}) - | awk -v n=${target_number} '$1<=n' > target.txt

#bash ${pipeline}/supp/edit-and-lift-input-annotations-to-assembly-annotations.sh ${target_assembly} ${target_annotations} | awk -f ${pipeline}/lift/lift-assembly-annotations-to-input-annotations.awk <(awk '$0~/^>/{len[$2]=$3;len[-$2]=$3;next}{counter++; for(i=1;i<=NF;i++){clen[counter]+=len[$i]}}END{for(i=1;i<=counter;i++){print ">"i, i, clen[i]};for(i=1;i<=counter;i++){print i}}' ${target_assembly}) - | awk -v n=${target_number} '$1<=n' > target.txt

#query
#bash ${pipeline}/supp/edit-and-lift-input-annotations-to-assembly-annotations.sh -b ${bundle_size} ${query_assembly} ${query_annotations} | awk -f ${pipeline}/lift/lift-assembly-annotations-to-input-annotations.awk <(awk '$0~/^>/{len[$2]=$3;len[-$2]=$3;next}{counter++; for(i=1;i<=NF;i++){clen[counter]+=len[$i]}}END{for(i=1;i<=counter;i++){print ">"i, i, clen[i]};for(i=1;i<=counter;i++){print i}}' ${query_assembly}) - | awk -v n=${query_number} '$1<=n' > query.txt

bash ${pipeline}/supp/edit-and-lift-input-annotations-to-assembly-annotations.sh ${query_assembly} ${query_annotations} | awk -f ${pipeline}/lift/lift-assembly-annotations-to-input-annotations.awk <(awk '$0~/^>/{len[$2]=$3;len[-$2]=$3;next}{counter++; for(i=1;i<=NF;i++){clen[counter]+=len[$i]}}END{for(i=1;i<=counter;i++){print ">"i, i, clen[i]};for(i=1;i<=counter;i++){print i}}' ${query_assembly}) - | awk -v n=${query_number} '$1<=n' > query.txt

# merge annotations into dotplot file && supplement with lengths
awk 'BEGIN{FS="\t"; OFS="\t"}FILENAME==ARGV[1]{counter[substr($8,2)]++; targetchr[substr($8,2)]=$1; targetpos[substr($8,2)]=$9; orientation[substr($8,2)]=substr($8,1,1); next}{counter[substr($8,2)]++; querychr[substr($8,2)]=$1; querypos[substr($8,2)]=$9; if (substr($8,1,1)==orientation[substr($8,2)]){orientation[substr($8,2)]="+"}else{orientation[substr($8,2)]="-"}}END{for(i in targetpos){if((i in querypos)&&counter[i]==2){n=split(i,a,":"); score=a[n]; print targetchr[i], targetpos[i], querychr[i], querypos[i], orientation[i], score}}}' <(awk '$0!="chr1	x1	x2	chr2	y1	y2	color	id	X1	X2	Y1	Y2"' target.txt) <(awk '$0!="chr1	x1	x2	chr2	y1	y2	color	id	X1	X2	Y1	Y2"' query.txt) | awk 'FILENAME==ARGV[1]{if($0~/^>/){tlen[$2]=$3;tlen[-$2]=$3; next}; tcounter++; for(i=1;i<=NF;i++){tchromlen[tcounter]+=tlen[$i]};next}FILENAME==ARGV[2]{if($0~/^>/){qlen[$2]=$3;qlen[-$2]=$3; next}; qcounter++; for(i=1;i<=NF;i++){qchromlen[qcounter]+=qlen[$i]};next}{FS="\t"; OFS="\t"}{print $0, tchromlen[$1], qchromlen[$3]}' ${target_assembly} ${query_assembly} - > target_vs_query.dotplot.txt

# cleanup
rm target.txt query.txt target_annotations.txt query_annotations.txt
