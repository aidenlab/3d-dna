#/bin/bash

pipeline=`cd "$( dirname $0)" && cd .. && pwd`

entry=$1

echo $entry

[ -d ${entry} ] || mkdir ${entry}

cd ${entry}

released="DNAZoo/50_mammals"
unpublished="DNAZoo/Unpublished"

dropbox_uploader.sh list ${released}"/"${entry}
if [ $? -eq 0 ]; then
    current=${released}
else
    dropbox_uploader.sh list ${unpublished}"/"${entry}
    [ $? -eq 0 ] && current=${unpublished}
fi

[ -z $current ] && { echo ":( Folder not found. Exiting!" && exit 1; }

dropbox_uploader.sh download ${current}"/"${entry}"/README.json" "README.json"

organism=$(jq '.organism | .binomial' README.json | awk '{gsub("\"","")}1')

echo "organism: "$organism

IFS='\n' read -r -a array <<<$(cat <(jq '.dnaSamples[] | .id' README.json) <(jq '.hicSamples[] | .id' README.json) | sort -u | awk '{gsub("\"","")}1')

if [ ${#array[@]} -eq 1 ];
then
    isolate=${array[0]}
else
    echo ":| Can't handle this case yet. Exiting!"
    exit 1
fi

echo "isolate: "$isolate

chromosomes=$(jq '.chromlengthAssembly | .karyotype' README.json | awk '{gsub("\"","");split($0,a,"="); if(a[1]=="2n"){print a[2]/2}else if (a[1]=="n"){print a[2]}}')

echo "chromosomes: "$chromosomes

if [ -z $chromosomes ] || [ "$chromosomes" == "" ]; then
    echo "No chromosome count listed. Exiting!"
    exit 1
fi

fasta=$(jq '.chromlengthAssembly | .name' README.json | awk '{gsub("\"","")}1')

dropbox_uploader.sh download ${current}"/"${entry}"/"$fasta".fasta.gz" $fasta".fasta.gz"

[ -f $fasta".fasta.gz" ] || { ":( Something went wrong. Exiting!" && exit 1; }

pigz -d $fasta".fasta.gz"

awk -v organism="${organism}" -v isolate="${isolate}" -v chromosomes="${chromosomes}" 'FILENAME==ARGV[1]{skip[$1]=1;next}$0~/^>/{test=1}($1 in skip){test=0}!test{next}$0!~/^>/{print; next}{counter++}{$0=$1" [organism="organism"] [isolate="isolate"]"}counter<=chromosomes{$0=$0" [location=chromosome] [chromosome="counter"]"}1' <(awk -v ignore_sorting=1 -f ${pipeline}/utils/generate-assembly-file-from-fasta.awk ${fasta}".fasta" | awk '$0~/^>/&&$3<200{print $1}') ${fasta}".fasta" > ${fasta}".mod_headers.fasta"

mv ${fasta}".mod_headers.fasta" ${fasta}".fasta"

pigz ${fasta}".fasta"