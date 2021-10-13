#!/bin/bash
## Wrapper script to generate final fasta as well as various component fastas (unprompted)
## Adds 500bp gaps between assembly components scaffolded via Hi-C
## TODO: make gap length a parameter
## Written by OD

pipeline=`cd "$( dirname $0)" && cd .. && pwd`

## unprompted
gap_size=500
subtiny_size=1000

# default
label="HiC"

USAGE="
***********************************************
./finalize-output.sh -c <number_of_chromosomes> -s <tiny_threshold> -g <gap_size> -l <label> -o <organism> -i <isolate> <cprops> <asm> <fasta> <type>
***********************************************
"

## HANDLE OPTIONS

while getopts "c:s:l:g:o:i:h" opt; do
case $opt in
	h) echo "$USAGE"
		exit 0
	;;
	c) re='^[0-9]+$'
		if [[ $OPTARG =~ $re ]] && [[ $OPTARG -gt 0 ]]; then
			echo "... -c flag was triggered, will add [location=chromosome] and [chromosome=<chrom_number>] to output sequence names." >&1
	        chrom_num="$OPTARG"
        else
	          echo ":( Wrong syntax for chromosome number. Exiting!" >&2 && exit 1
  	  	fi
	;;
	s) 	re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]] && [[ $OPTARG -gt 0 ]]; then
	          echo "... -s flag was triggered, treating all contigs/scaffolds shorter than $OPTARG as unattempted." >&1
	          input_size="$OPTARG"
        else
	          echo ":( Wrong syntax for minimal input contig/scaffold size. Exiting!" >&2 && exit 1
  	  	fi
	;;
	g) 	re='^[0-9]+$'
		if [[ $OPTARG =~ $re ]] && [[ $OPTARG -gt 0 ]]; then
			echo "... -g flag was triggered, making gap size between scaffolded draft sequences to be equal to $OPTARG." >&1
			gap_size="$OPTARG"
		else
			echo ":( Wrong syntax for default gap size parameter value. Using default gap_size=${gap_size}!" >&2
	fi
	;;
	o)
		echo "... -o flag was triggered, will add [organism="$OPTARG"] to output sequence names." >&1
		organism="$OPTARG"
	;;
	i)
		echo "... -i flag was triggered, will add [isolate="$OPTARG"] to output sequence names." >&1
		isolate="$OPTARG"
	;;
	l) 
		echo "... -l flag was triggered. Output sequence will will appear as ${OPTARG}_HiC.fasta."
		label="$OPTARG"
	;;
	*) echo "$USAGE" >&2
		exit 1
	;;
esac
done
shift $(( OPTIND-1 ))

## HANDLE ARGUMENTS

[ -z $1 ] || [ -z $2 ] || [ -z $3 ] || [ -z $4 ] && echo >&2 ":( Some input seems to be missing." && echo >&2 "$USAGE" && exit 1

cprops=$1
asm=$2
fasta=$3
type=$4

case $type in
	"draft")
		[ -z $input_size ] && echo >&2 ":( Please type in the minimal contig/scaffold input size -s <num>. Exiting!" && exit 1
		echo "Analyzing the draft assembly"
		
		# Note that if the draft is not trimmed might need to trim first...
		
		awk -v except=1 -f ${pipeline}/finalize/make-fasta-subset.awk <( awk -v subtiny=${subtiny_size} '$3<subtiny' ${cprops}) ${fasta} > draft.fasta
		awk -v except=1 -f ${pipeline}/finalize/make-fasta-subset.awk <( awk -v tiny=${input_size} '$3<tiny' ${cprops}) ${fasta} > input.fasta
		awk -f ${pipeline}/finalize/make-fasta-subset.awk <( awk -v subtiny=${subtiny_size} -v tiny=${input_size} '$3>=subtiny&&$3<tiny' ${cprops}) ${fasta} > tiny.fasta
		ln -sf tiny.fasta unattempted.fasta		
	;;
	"raw")
		echo "Analyzing the tiled assembly"
		[ -z ${input_size} ] || [ -z ${chrom_num} ] && echo >&2 ":( Please type in the number of chromosomes -c <num> and the tiny size -s <num>. Exiting!" && exit 1
		
		awk 'gsub("{||}","")' ${asm} > temp.asm
		asm=temp.asm
		bash ${pipeline}/finalize/remove-N-overhangs-from-asm.sh ${cprops} ${asm} ${fasta}
		
		prefix=`basename ${cprops} .cprops`
		cprops=${prefix}.no_overhangs.cprops
		prefix=`basename ${asm} .asm`
		asm=${prefix}.no_overhangs.asm
		prefix=`basename ${fasta} .fa`
		prefix=`basename ${prefix} .fna`
		prefix=`basename ${prefix} .fasta`
		fasta=${prefix}.no_overhangs.fasta
		
		awk -v chrom_num=${chrom_num} 'FILENAME==ARGV[1]{oname[$2]=$1;next}FNR<=chrom_num{gsub("-","");for(i=1;i<=NF;i++){print oname[$i]}}' ${cprops} ${asm} | awk -f ${pipeline}/finalize/make-fasta-subset.awk - ${fasta} > resolved.fasta
		awk -v chrom_num=${chrom_num} -v tiny=${input_size} 'FILENAME==ARGV[1]{oname[$2]=$1;split($1,b,":::overhang_||:::gap");len[b[1]]+=$3;next}FNR>chrom_num{gsub("-","");for(i=1;i<=NF;i++){split(oname[$i],b,":::overhang_||:::gap");if(oname[$i]~/:::fragment_/||len[b[1]]>=tiny){print oname[$i]}}}' ${cprops} ${asm} | awk -f ${pipeline}/finalize/make-fasta-subset.awk - ${fasta} > unresolved-and-inconsistent.fasta
		
		rm temp.asm
				
	;;
	"final")
	
		echo "Analyzing the assembly"

		# trim N overhangs
		echo "...trimming N overhangs"
	
		bash ${pipeline}/finalize/remove-N-overhangs-from-asm.sh ${cprops} ${asm} ${fasta}
		prefix=`basename ${cprops} .cprops`
		cprops=${prefix}.no_overhangs.cprops
		prefix=`basename ${asm} .asm`
		asm=${prefix}.no_overhangs.asm
		prefix=`basename ${fasta} .fa`
		prefix=`basename ${prefix} .fna`
		prefix=`basename ${prefix} .fasta`
		fasta=${prefix}.no_overhangs.fasta

		# riffle
		echo "...adding gaps"
		
		gap_id=`awk 'END{print \$2+1}' ${cprops}`
	
		awk -v riffle=${gap_id} -f ${pipeline}/finalize/riffle-asm.awk ${asm} > temp.asm
	
		cp ${cprops} temp.cprops
		echo "hic_gap_${gap_id} ${gap_id} ${gap_size}" >> temp.cprops
	
		cp ${fasta} temp.fasta
		echo ">hic_gap_${gap_id}" >> temp.fasta
		awk -v gap_size=${gap_size} 'BEGIN{for(i=1; i<=gap_size;i++){str=str"N"}; print str}' >> temp.fasta

		bash ${pipeline}/finalize/construct-fasta-from-asm.sh temp.cprops temp.asm temp.fasta | awk -f ${pipeline}/utils/wrap-fasta-sequence.awk - | awk -v organism="${organism}" -v isolate="${isolate}" -v chrom_num=${chrom_num} '$0!~/>/{print; next}organism{$1=$1" [organism="organism"]"}isolate{$1=$1" [isolate="isolate"]"}{counter++}(counter<=chrom_num){$1=$1" [location=chromosome] [chromosome="counter"]"}1' > ${label}_HiC.fasta
		
		# clean up: remove no_overhangs files
		
		awk -v disable_checks=1 -f ${pipeline}/utils/convert-cprops-and-asm-to-assembly.awk temp.cprops temp.asm && mv temp.assembly ${label}_HiC.assembly
		rm temp.cprops temp.asm temp.fasta ${cprops} ${asm} ${fasta}
		
		exit 
		
		# remove rest from script, not specialized script for this
# skip component analysis if no optional data like chrom number
		[ -z ${input_size} ] || [ -z ${chrom_num} ] || [ -z ${label} ] && echo >&2 ":| Do not know the number of expected chromosomes: skipping component fasta analysis. Done!" && rm temp.cprops temp.asm temp.fasta && exit 0

		awk -v chrom_num=${chrom_num} 'NR<=chrom_num' temp.asm > temp.chr-length.asm
	
		awk -v subtiny=${subtiny_size} 'FILENAME==ARGV[1]{split($1,a,":::overhang_"); len[a[1]]+=$3; oname[$2]=a[1]; next}{gsub("-","")}(NF==1&&oname[$1]!~/:::fragment_/&&len[oname[$1]]<subtiny){print}' temp.cprops temp.asm > temp.subtiny.asm
	
		awk -v subtiny=${subtiny_size} -v tiny=${input_size} 'FILENAME==ARGV[1]{split($1,a,":::overhang_"); len[a[1]]+=$3; oname[$2]=a[1]; next}{gsub("-","")}(NF==1&&oname[$1]!~/:::fragment_/&&len[oname[$1]]>=subtiny&&len[oname[$1]]<tiny){print}' temp.cprops temp.asm > temp.tiny.asm
	
		cat temp.chr-length.asm temp.subtiny.asm temp.tiny.asm | awk 'FILENAME==ARGV[1]{skip[$0]=1;next}(!skip[$0]){print}' - temp.asm > temp.small.asm

		# make scaffold sequences of groups and calculate statistics
		echo "...generating fastas"
		

		bash ${pipeline}/finalize/construct-fasta-from-asm.sh temp.cprops temp.chr-length.asm temp.fasta | awk -f ${pipeline}/utils/wrap-fasta-sequence.awk - > chr-length.fasta
				
		awk '$0~/>/{if(c){print c}; c=0; next}{c+=length}END{print c}' chr-length.fasta > T3.dat
		awk '$0~/>/{if(c){print c}; c=0; next}{gsub("N","");gsub("n","");c+=length}END{print c}' chr-length.fasta > tmp
#		paste T3.dat tmp > T3.dat.tmp && mv T3.dat.tmp T3.dat && rm tmp
		
		bash ${pipeline}/finalize/construct-fasta-from-asm.sh temp.cprops temp.small.asm  temp.fasta | awk -f ${pipeline}/utils/wrap-fasta-sequence.awk - > small.fasta
		
		bash ${pipeline}/finalize/construct-fasta-from-asm.sh temp.cprops temp.tiny.asm  temp.fasta | awk -f ${pipeline}/utils/wrap-fasta-sequence.awk - > tiny.fasta
	
		# merge final output
		
		cat chr-length.fasta small.fasta | awk -v label=${label} '$0~/>/{counter++; $0=">"label"_hic_scaffold_"counter}1' > ${label}_HiC.from_input.fasta
		
		cat chr-length.fasta small.fasta tiny.fasta | awk -v label=${label} '$0~/>/{counter++; $0=">"label"_hic_scaffold_"counter}1' > ${label}_HiC.from_draft.fasta
		
		cat small.fasta tiny.fasta > small-and-tiny.fasta
	
		# clean up: remove no_overhangs files
		rm ${cprops} ${asm} ${fasta}
		# clean up: remove component cprops and asm
		rm temp.fasta temp.cprops temp.asm temp.small.asm temp.tiny.asm temp.chr-length.asm temp.subtiny.asm
		# clean up: comment if component fastas are needed
		rm chr-length.fasta small.fasta tiny.fasta small-and-tiny.fasta

	;;
	*)
		echo >&2 ":( Unknown type. Please choose one of the following: draft/tiled/merged. Exiting!" && exit 1
	;;
esac
