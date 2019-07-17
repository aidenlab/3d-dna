#!/bin/bash

#### Description: Wrapper script to polish any given assembly. Polish is a process in which a mismatch detector is run and pieces of an assembly between mismatches are treated as an input. Polish cannot insert pieces unlike a typical iterative assembly step but it can avoid some pitfalls caused by very small or unusually behaving contigs and scaffolds. The idea is that you typically run mismatch detector at a pretty large scale to make sure there are no gaping global misassemblies left. To control what is going on one might want to check if the number of mismatches pre to post-polish has decreased (currently disabled).
#### Usage: run-asm-polisher.sh -j <path_to_current_hic_file> -a <path_to_scaf_annotation_file> -b <path_to_superscaf_annotation_file> -w <coarse_res_for_mismatch> -n <fine_res_for_mismatch> -d <depletion_region> <path_to_original_cprops> <path_to_original_mnd_file> <path_to_current_cprops> <path_to_current_asm_file>
#### Input: cprops and mnd for original input contigs/scaffolds; cprops and asm for current input contigs/scaffolds.
#### Optional input: Juicebox .hic file for current assembly (-j option).
#### Output: cprops and asm of the polished assembly. Additional files include the new mnd and the new .hic files.
#### Parameters: primarily those that will be passed to the mismatch detector.
#### Unprompted: -p for use of GNU Parallel; -c for percent to saturate, -k for sensitivity, -b for balancing type (mismatch detector unprompted parameters).
#### Dependencies: mismatch-detector, editor, scaffolder, polish-specific files (edit-asm-according-to-new-cprops.sh).
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu. Version dated 01/22/2017

## Set defaults
input_size=1000000
wide_bin=100000
wide_depletion=3000000
narrow_bin=1000

## Set unprompted defaults
use_parallel=true		# use GNU Parallel to speed-up calculations (default)
mapq=1				# minimal mapping quality
k=55					# sensitivity to depletion score (50% of expected is labeled as a mismatch)
pct=5					# default percent of map to saturate
norm="KR"				# use an unbalanced contact matrix for analysis

## HANDLE OPTIONS

while getopts "hs:j:a:b:w:n:d:k:c:b:p:q:" opt; do
case $opt in
    h) echo "$USAGE" >&1
        exit 0
    ;;
    p)	if [ $OPTARG == true ] || [ $OPTARG == false ]; then
    	    echo ":) -p flag was triggered. Running with GNU Parallel support parameter set to $OPTARG." >&1
			use_parallel=$OPTARG
    	else
    		echo ":( Unrecognized value for -p flag. Running with default parameters (-p true)." >&2
    	fi
    ;;
    q)  re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]]; then
            echo ":) -q flag was triggered, performing polishing taking into account signal with minimum $OPTARG mapping quality" >&1
            mapq=$OPTARG
        else
            echo ":( Wrong syntax for mapping quality size. Using the default value $mapq" >&2
        fi
    ;;
    s)  re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]]; then
            echo "...-s flag was triggered, will ignore all scaffolds shorter than $OPTARG for polishing" >&1
            input_size=$OPTARG
        else
            echo ":( Wrong syntax for input size. Using the default value ${input_size}" >&2
        fi
    ;;
   	j) if [ -s $OPTARG ]; then 
			echo "...-j flag was triggered, will use Juicebox map $OPTARG" >&1
			current_hic=$OPTARG
		else
			echo ":( Juicebox file not found. Will run visualize script from scratch" >&2
		fi	
	;;
   	a) if [ -s $OPTARG ]; then 
			echo "...-a flag was triggered, will use scaffold annotation file $OPTARG" >&1
			current_scaf=$OPTARG
		else
			echo ":( Scaffold annotation file not found. Will run visualize script from scratch" >&2
		fi	
	;;
	b) if [ -s $OPTARG ]; then 
			echo "...-b flag was triggered, will use superscaffold annotation file $OPTARG" >&1
			current_superscaf=$OPTARG
		else
			echo ":( Superscaffold annotation file not found. Will run visualize script from scratch" >&2
		fi	
	;;
    w)  re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]]; then
            echo ":) -w flag was triggered, performing cursory search for mismatches at $OPTARG resolution" >&1
            wide_bin=$OPTARG
        else
            echo ":( Wrong syntax for bin size. Using the default value 25000" >&2
        fi
    ;;
    n)  re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]]; then
            echo ":) -n flag was triggered, performing mismatch region thinning at $OPTARG resolution" >&1
            narrow_bin=$OPTARG
        else
            echo ":( Wrong syntax for mismatch localization resolution. Using the default value 1000" >&2
        fi
    ;;
	d)  re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]]; then
            echo ":) -d flag was triggered, depletion score will be averaged across a region bounded by $OPTARG superdiagonal" >&1
            wide_depletion=$OPTARG
        else
            echo ":( Wrong syntax for mapping quality. Using the default value dep_size=100000" >&2
        fi
    ;;
    k)  re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]] && [[ $OPTARG -gt 0 ]] && [[ $OPTARG -lt 100 ]]; then
            echo ":) -k flag was triggered, starting calculations with ${OPTARG}% depletion as mismatch threshold" >&1
            k=$OPTARG
        else
            echo ":( Wrong syntax for mismatch threshold. Using the default value k=50" >&2
        fi
    ;;
    c)  re='^[0-9]+\.?[0-9]*$'
if [[ $OPTARG =~ $re ]] && [[ ${OPTARG%.*} -ge 0 ]] && ! [[ "$OPTARG" =~ ^0*(\.)?0*$ ]] && [[ $((${OPTARG%.*} + 1)) -le 100 ]]; then
        	echo ":) -c flag was triggered, starting calculations with ${OPTARG}% saturation level" >&1
        	pct=$OPTARG
        else
        	echo ":( Wrong syntax for saturation threshold. Using the default value pct=${pct}" >&2
        fi
	;;
    b)	if [ $OPTARG == NONE ] || [ $OPTARG == VC ] || [ $OPTARG == VC_SQRT ] || [ $OPTARG == KR ]; then
    	    echo ":) -b flag was triggered. Type of norm chosen for the contact matrix is $OPTARG." >&1
			norm=$OPTARG
    	else
    		echo ":( Unrecognized value for -b flag. Running with default parameters (-b NONE)." >&2
    	fi
    ;;
    *) echo "$USAGE" >&2
        exit 1
    ;;
esac
done

shift $(( OPTIND-1 ))

## HANDLE ARGUMENTS: TODO check file format
if [ $# -lt 4 ]; then
    echo ":( Required arguments not found. Please double-check your input!!" >&2
    echo "$USAGE" >&2
    exit 1
fi

orig_cprops=$1
orig_mnd=$2
current_cprops=$3
current_asm=$4

if [ ! -f ${orig_cprops} ] || [ ! -f ${orig_mnd} ] || [ ! -f ${current_cprops} ] || [ ! -f ${current_asm} ] ; then
	echo >&2 ":( Required files not found. Please double-check your input!!" && exit 1
fi

## CHECK DEPENDENCIES

pipeline=`cd "$( dirname $0)" && cd .. && pwd`

## PREP

id=`basename ${current_cprops} .cprops`

STEP="polish"

# PREPARATORY: split asm into resolved and unresolved pieces:
awk 'NR==1' ${current_asm} > temp_resolved.asm
awk 'NR>1' ${current_asm} > temp_unresolved.asm

## MAIN FUNCTION

#	0) Check if Juicebox file for the current assembly has been passed. If not build a map. TODO: make sure that the required resolution is included. TODO: Enable at some point.

if [ -z ${current_hic} ] || [ -z ${current_scaf} ] || [ -z ${current_superscaf} ] ; then
	echo "...no hic file and/or annotation files have been provided with input, building the hic map from scratch" >&1
	## TODO: check w/o parallel
	bash ${pipeline}/edit/edit-mnd-according-to-new-cprops.sh ${current_cprops} ${orig_mnd} > `basename ${current_cprops} .cprops`.mnd.txt
	current_mnd=`basename ${current_cprops} .cprops`.mnd.txt
	bash ${pipeline}/visualize/run-asm-visualizer.sh -q ${mapq} -i -c ${current_cprops} ${current_asm} ${current_mnd}
	current_hic=`basename ${current_asm} .asm`.hic
	current_scaf=`basename ${current_asm} .asm`_asm.scaffold_track.txt
	current_superscaf=`basename ${current_asm} .asm`_asm.superscaf_track.txt
fi

#	1) Annotate mismatches in current assembly
bash ${pipeline}/edit/run-mismatch-detector.sh -p ${use_parallel} -c ${pct} -w ${wide_bin} -k ${k} -d ${wide_depletion} -n ${narrow_bin} ${current_hic}

# store intermediate mismatch stuff	- not necessary
mv depletion_score_wide.wig ${id}.${STEP}.depletion_score_wide.wig
mv depletion_score_narrow.wig ${id}.${STEP}.depletion_score_narrow.wig
mv mismatch_wide.bed ${id}.${STEP}.mismatch_wide.bed
mv mismatch_narrow.bed ${id}.${STEP}.mismatch_narrow.bed

# convert bed track into 2D annotations
resolved=$(awk 'NR==2{print $3}' ${current_superscaf})	#scaled coordinates
awk -v bin_size=${narrow_bin} -f ${pipeline}/edit/overlay-edits.awk ${current_scaf} ${id}.${STEP}.mismatch_narrow.bed | awk -v r=${resolved} 'NR==1||$3<=r' > ${id}.${STEP}.suspect_2D.txt

# split into mismatches and edits
	awk 'NR==1||$8=="mismatch"' ${id}.${STEP}.suspect_2D.txt > ${id}.${STEP}.mismatches_2D.txt
	awk 'NR==1||$8=="debris"' ${id}.${STEP}.suspect_2D.txt > ${id}.${STEP}.edits_2D.txt

#	2) Deal with mismatches: break resolved asm at joints associated with mismatches	
	awk -f ${pipeline}/lift/lift-asm-annotations-to-input-annotations.awk ${current_cprops} ${current_asm} ${id}.${STEP}.mismatches_2D.txt | awk 'FILENAME==ARGV[1]{cname[$1]=$2;next}FNR>1{if($2==0){type[cname[$1]]="h"}else{type[cname[$1]]="t"}}END{for(i in type){print i, type[i]}}' ${current_cprops} - | awk 'FILENAME==ARGV[1]{type[$1]=$2; if($2=="h"){type[-$1]="t"}else{type[-$1]="h"}; next}{str=""; for(i=1;i<=NF;i++){if($i in type){if (type[$i]=="h"){print str; str=$i}else{str=str" "$i; print str; str=""}}else{str=str" "$i}}; print str}' - temp_resolved.asm | sed '/^$/d' | awk '{$1=$1}1' > temp_resolved.asm.new && mv temp_resolved.asm.new temp_resolved.asm
		
#	3) Apply edits
	
	# reconstruct the edits file

	awk 'BEGIN{OFS="\t"; print "chr1", "x1", "x2", "chr2", "y1", "y2", "color", "id", "X1", "X2", "Y1", "Y2"}$1~/:::debris/{print $1, 0, $3, $1, 0, $3, "0,0,0", "debris", 0, $3, 0, $3}' ${current_cprops} | awk -f ${pipeline}/lift/lift-input-annotations-to-asm-annotations.awk ${current_cprops} <(awk '{print $2}' ${current_cprops}) - | awk -f ${pipeline}/lift/lift-asm-annotations-to-input-annotations.awk ${orig_cprops} <(awk '{print $2}' ${orig_cprops}) - > temp.pre_polish_edits.txt
	
	current_edits=temp.pre_polish_edits.txt
	bash ${pipeline}/lift/lift-edit-asm-annotations-to-original-input-annotations.sh ${orig_cprops} ${current_cprops} ${current_asm} ${id}.${STEP}.edits_2D.txt > h.edits.txt
	awk 'NR==1' "h.edits.txt" > temp
	{ awk 'NR>1' ${current_edits} ; awk 'NR>1' "h.edits.txt" ; } | sort -k 1,1 -k 2,2n >> temp
	mv temp temp.post_polish_edits.txt
	polish_edits=temp.post_polish_edits.txt
	bash ${pipeline}/edit/apply-edits-prep-for-next-round.sh -p ${use_parallel} -r ${STEP} ${polish_edits} ${orig_cprops} ${orig_mnd}
	mv `basename ${orig_cprops} .cprops`.${STEP}.cprops $id.${STEP}.cprops
	mv `basename ${orig_mnd} .txt`.${STEP}.txt $id.${STEP}.mnd.txt
	polish_cprops=$id.${STEP}.cprops
	polish_mnd=$id.${STEP}.mnd.txt
	
#	4) Lift current assembly to new cprops
	bash ${pipeline}/edit/edit-asm-according-to-new-cprops.sh ${polish_cprops} ${current_cprops} temp_resolved.asm > new_resolved.asm
	bash ${pipeline}/edit/edit-asm-according-to-new-cprops.sh ${polish_cprops} ${current_cprops} temp_unresolved.asm > new_unresolved.asm
	
#	5) Prepare for polish run: break resolved at debris, filter pieces smaller than input_size
	awk -v input_size=${input_size} 'function printout(str){if(c>=input_size){print substr(str,2)>"h.scaffolds.original.notation.step.0.txt"}else{print substr(str,2)>"h.dropouts.step.0.txt"}}FILENAME==ARGV[1]{len[$2]=$3; len[-$2]=$3; if($1~/:::debris/){remove[$2]=1; remove[-$2]=1}; next}{str=""; for(i=1;i<=NF;i++){if($i in remove){if(str!=""){printout(str)}; print $i > "h.dropouts.step.0.txt"; str=""; c=0}else{str=str" "$i; c+=len[$i]}}; if(str!=""){printout(str)}}' ${polish_cprops} new_resolved.asm
	cat new_unresolved.asm >> h.dropouts.step.0.txt

mv h.dropouts.step.0.txt do_not_delete.dropouts.step.0.txt

#	6) Run TIGer
	bash ${pipeline}/scaffold/run-tiger-scaffolder.sh -p ${use_parallel} -q ${mapq} -s ${input_size} ${polish_cprops} ${polish_mnd}
	polish_asm=`basename ${polish_cprops} .cprops`.asm

mv do_not_delete.dropouts.step.0.txt h.dropouts.step.0.txt
mv ${polish_asm} h.scaffolds.original.notation.step.0.txt

#	6) Run LIGer (for things TIGer was not able to join - not necessary, but for megascaffold consistency)
	bash ${pipeline}/scaffold/run-liger-scaffolder.sh -p ${use_parallel} -q ${mapq} -s ${input_size} ${polish_cprops} ${polish_mnd}
	polish_asm=`basename ${polish_cprops} .cprops`.asm

#	7) Visualize output
	bash ${pipeline}/visualize/run-asm-visualizer.sh -p ${use_parallel} -q ${mapq} -i -c ${polish_cprops} ${polish_asm} ${polish_mnd}
	
#	8) Cleanup
	rm ${polish_mnd} temp_resolved.asm temp_unresolved.asm temp.pre_polish_edits.txt temp.post_polish_edits.txt new_resolved.asm new_unresolved.asm
