#!/bin/bash

USAGE="
	./snp-mnd-to-diploid-mnd.sh <phased_vcf> <snp_mnd_file>

"

vcf=$1
mnd=$2

pipeline=`cd "$( dirname $0)" && cd .. && pwd`

#awk '$1!=prev&&NR>1&&$NF~/|/{max=0; for(i in c){if(c[i]>max){max=c[i]; maxi=i}}; print prev" "maxi; delete c;}{prev=$1; gsub(/.+:/,"",$NF); c[$NF]++}END{max=0; for(i in c){if(c[i]>max){max=c[i]; maxi=i}}; print prev" "maxi;}' $vcf | awk 'FILENAME==ARGV[1]{split($0,a," "); largestblock[a[1]":0|1:"a[2]]=1; largestblock[a[1]":1|0:"a[2]]=1; next}FILENAME==ARGV[2]{if(!($1":"$NF in largestblock)){next}; if($NF~/0\|1/){allele[$1":"$2":"$4]=$1"r"; allele[$1":"$2":"$5]=$1"a"}else{allele[$1":"$2":"$4]=$1"a"; allele[$1":"$2":"$5]=$1"r"}; pos[$1":"$2":"$4]=$2; pos[$1":"$2":"$5]=$2; next}(!($NF in dup)){if(allele[$2]<=allele[$6]){$3=pos[$2]; $2=allele[$2]; $7=pos[$6]; $6=allele[$6]}else{$3=pos[$6]; $7=pos[$2]; tmp=$2; $2=allele[$6]; $6=allele[tmp]; tmp=$1; $1=$5; $5=tmp; tmp=$4; $4=$8; $8=tmp};if($2!=""&&$3!=""&&$6!=""&&$7!=""){print $1, $2, $3, $4, $5, $6, $7, $8; dup[$NF]=1}}' - $vcf $mnd | sort -k2,2d -k6,6d --parallel=48 -S32G > diploid.mnd.txt


awk '$1!=prev&&NR>1&&$NF~/|/{max=0; for(i in c){if(c[i]>max){max=c[i]; maxi=i}}; print prev" "maxi; delete c;}{prev=$1; gsub(/.+:/,"",$NF); c[$NF]++}END{max=0; for(i in c){if(c[i]>max){max=c[i]; maxi=i}}; print prev" "maxi;}' $vcf | awk 'FILENAME==ARGV[1]{split($0,a," "); largestblock[a[1]":0|1:"a[2]]=1; largestblock[a[1]":1|0:"a[2]]=1; next}FILENAME==ARGV[2]{if(!($1":"$NF in largestblock)){next}; if($NF~/0\|1/){allele[$1":"$2":"$4]=$1"r"; allele[$1":"$2":"$5]=$1"a"}else{allele[$1":"$2":"$4]=$1"a"; allele[$1":"$2":"$5]=$1"r"}; next}{if(allele[$2]!=""){$2=allele[$2]; $6=$2}else{$6=allele[$6]; $2=$6}; if($2!=""){print $1, $2, $3, $4, $5, $6, $7, $8}}' - $vcf $mnd | sort -k2,2d -k6,6d --parallel=48 -S32G > diploid.mnd.txt