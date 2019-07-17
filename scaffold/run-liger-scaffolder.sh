#!/bin/bash
## LIGer scaffolder wrapper script
## Written by Olga Dudchenko and Sanjit Batra

USAGE="
*****************************************************
This is a wrapper for a Hi-C Limitless Iterative Greedy genome assembly (LIGer) algorithm, version date: Dec 7, 2016.

Usage: ./run-liger-scaffolder.sh [-h] [-s minimal_scaffold_size] [-t link_threshold] [-q mapq] path_to_cprops_file path_to_merge_nodups_file

ARGUMENTS:
path_to_cprops_file     Path to (prefiltered) cprops file listing contigs for which LIGer scaffolding will be attempted
path_to_merge_nodups_file  Path to merge_nodups Juicer output file

OPTIONS:
-h                      Shows this help
-s size					Set minimal contig/scaffold size to use as input
-q mapq                 Set threshold for Hi-C reads mapping quality (default is 1)
-p true/false			Use GNU Parallel to speed up calculations (default is true)
-t link_threshod      	Set threshold for joining links [not working yet, uses default]

Uses scrape-mnd.awk, generate-unsorted-confidence-table.awk, confidence-to-assembly.awk, scaffolds-to-original-notation.awk and drop-smallest-dubious-element.awk that should be in the same folder as the wrapper script.

In the current version ordering and orienting contigs at each iteration is based on the data from  whole inter-contig contact matrices. At each iteraction the contigs and corresponding matrices are updated, the idea being that some low-confidence scaffolds links may benefit from more data provided by the already-joined contigs and get resolved.

Note that in the current version the input is expected in the cprops format.
*****************************************************
"

## Set current defaults
SIZE=15000
MAPQ=1
use_parallel=true

## HANDLE OPTIONS
while getopts "s:q:p:t:h" opt; do
case $opt in
    h) echo "$USAGE" >&2
        exit 0
    ;;
    s)  re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]] && [[ $OPTARG -gt 0 ]]; then
            echo ":) -s flag was triggered, starting calculations with $OPTARG threshold starting contig/scaffold size" >&1
            SIZE=$OPTARG
        else
            echo ":( Wrong syntax for minimal input contig/scaffold size. Using the default value SIZE=$SIZE" >&2
        fi
    ;;
    q)  re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]]; then
            echo ":) -q flag was triggered, starting calculations with $OPTARG threshold mapping quality" >&1
            MAPQ=$OPTARG
        else
            echo ":( Wrong syntax for mapping quality. Using the default value MAPQ=$MAPQ" >&2
        fi
    ;;
    p)	if [ $OPTARG == true ] || [ $OPTARG == false ]; then
    	    echo ":) -p flag was triggered. Running LIGer with GNU Parallel support parameter set to $OPTARG." >&1
			use_parallel=$OPTARG
    	else
    		echo ":( Unrecognized value for -p flag. Running LIGer with default parameters (-p true)." >&2
    	fi
    ;;
    t) echo "...-thr flag was triggered. Sorry, option not functional yet: using default thresholds" >&2
    ;;
    *) echo "$USAGE" >&2
        exit 1
    ;;
esac
done

shift $(( OPTIND-1 ))

## HANDLE ARGUMENTS: TODO check file format
if [ $# -lt 2 ]; then
    echo ":( Required arguments not found. Please double-check your input!!" >&2
    echo "$USAGE" >&2
    exit 1
fi
# Handle arguments: cprops file
contigPropFile=$1
if [[ ! -s "$contigPropFile" ]]; then
    echo ":( Cprops file not found. Please double-check your input!!" >&2
    echo "$USAGE" >&2
    exit 1
else
    # TODO: check that cprops file is in proper format
    echo "...Using cprops file: $contigPropFile"
fi

# Handle arguments: merged_nodups files
if [ $# -eq 2 ]; then
    # TODO: check that merged_nodups file is in proper format
    mergelib=$2
    echo "...Using merged_nodups file: $mergelib"
fi


## CHECK DEPENDENCIES

if [ $use_parallel == true ]; then
	type parallel >/dev/null 2>&1 || { echo >&2 ":( GNU Parallel support is set to true (default) but GNU Parallel is not in the path. Please install GNU Parallel or set -p option to false. Exiting!"; exit 1; }
fi

path_to_scripts=`cd "$( dirname $0)" && pwd`

scrape_contacts_script="$path_to_scripts""/scrape-mnd.awk"
merge_scores_script="$path_to_scripts""/merge-scores.awk"
compute_confidences_script="$path_to_scripts""/generate-unsorted-confidence-table.awk"
accept_links_script="$path_to_scripts""/confidence-to-assembly.awk"
update_assembly_script="$path_to_scripts""/scaffolds-to-original-notation.awk"
drop_dubious_script="$path_to_scripts""/drop-smallest-dubious-element.awk"


if [ ! -f $scrape_contacts_script ] || [ ! -f $merge_scores_script ] || [ ! -f $compute_confidences_script ] || [ ! -f $accept_links_script ] || [ ! -f $update_assembly_script ] || [ ! -f $drop_dubious_script ]; then
    echo ":( Relevant dependency scripts not found in bin folder. Exiting!" >&2
    exit 1
fi

## PREP FOR FIRST STEP - TODO: make restarting from any step possible
if [ ! -f "h.scaffolds.original.notation.step.""0"".txt" ]; then
	echo "...Scaffolding all scaffolds and contigs greater or equal to $SIZE bp."	
# thinking of introducing an unattempted flag, would influence things here
	gawk -v SIZE=${SIZE} '$3>=SIZE && $1!~/:::debris$/{print $2; next}{print $2 >"/dev/stderr"}' $contigPropFile > "h.scaffolds.original.notation.step.""0"".txt" 2>"h.dropouts.step.""0"".txt"
else
	echo "...Explicit scaffold set has been listed as input. Using set as a first iteration."
fi

STEP=1
echo "...Starting iteration # $STEP"

# MAIN LOOP

while true; do

	#do not enter the next iteration if nothing to assemble
	if [[ $(wc -l <"h.scaffolds.original.notation.step.""$((STEP-1))"".txt") -eq 1 ]]; then
		STEP=$((-1+$STEP))
		break
	fi

    #extract, relable and count reads from merged-nodups [TODO: rethink this part once mnd is deprecated]
	if [ $use_parallel == true ]; then
		parallel -a $mergelib --will-cite --jobs 80% --pipepart --block 1G "gawk -v MAPQ=$MAPQ -f $scrape_contacts_script $contigPropFile h.scaffolds.original.notation.step.$(($STEP-1)).txt - " | LC_ALL=C sort -k1,1 -k2,2 -k3,3n -s | gawk -f ${merge_scores_script} $contigPropFile "h.scaffolds.original.notation.step.""$(($STEP-1))"".txt" - > "h.scores.step.""$STEP"".txt"
	else
		gawk -v MAPQ="$MAPQ" -f $scrape_contacts_script $contigPropFile "h.scaffolds.original.notation.step.""$(($STEP-1))"".txt" "$mergelib" | gawk -f ${merge_scores_script} $contigPropFile "h.scaffolds.original.notation.step.""$(($STEP-1))"".txt" - > "h.scores.step.""$STEP"".txt"
	fi
	    
    #consolidate scrape data into double-sorted-confidence file
    gawk -f $compute_confidences_script "h.scores.step.""$STEP"".txt" | sort -r -gk4 -gk5 -S8G --parallel=48 -s > "h.double.sorted.confidence.step.""$STEP"".txt"

    #create new links between contigs based on confidence file
    gawk -f $accept_links_script "h.double.sorted.confidence.step.""$STEP"".txt" > "h.scaffolds.step.""$STEP"".txt"

    #update assembly file given the new links set
    gawk -f $update_assembly_script "h.scaffolds.step.""$STEP"".txt" "h.scaffolds.original.notation.step.""$(($STEP-1))"".txt" > "h.scaffolds.original.notation.step.""$STEP"".txt"

    #move to next step if the assembly was updated
    if [[ -s "h.scaffolds.step.""$STEP"".txt" ]]; then
        STEP=$((1+$STEP))
        echo "...Starting iteration # $STEP"

    #handle case when there are no links to add
    else
        ## TODO: is this the ultimate end?

        #try to push further by dropping dubious contigs

        if [ -f "h.dropouts.step.""$STEP"".txt" ]; then
            rm -f "h.dropouts.step.""$STEP"".txt"
        fi

        while true; do
            # Choose the smallest dubious contig to be dropped
            drop=$(gawk '($4==1){a[$1]; a[$2]}END{for (tmp in a) print tmp}' "h.double.sorted.confidence.step.""$STEP"".txt" | gawk -f $drop_dubious_script - $contigPropFile "h.scaffolds.original.notation.step.""$STEP"".txt")

            # One of ultimate end scenarios, probably there are better ways to get out
            if [ -z "$drop" ]; then
                break 2
            fi
            # split and overwrite input
            sed -n "$drop""p" "h.scaffolds.original.notation.step.""$STEP"".txt" >> "h.dropouts.step.""$STEP"".txt"
            gawk -v DROP="$drop" 'NR!=DROP' "h.scaffolds.original.notation.step.""$STEP"".txt" > "h.scaffolds.original.notation.tmp.step.""$STEP"".txt"
            mv "h.scaffolds.original.notation.tmp.step.""$STEP"".txt" "h.scaffolds.original.notation.step.""$STEP"".txt"

            # fix and overwrite scores
            gawk -v DROP="$drop" '($1!=DROP)&&($2!=DROP){if ($1>DROP) $1--; if ($2>DROP) $2--; print}' "h.scores.step.""$STEP"".txt" > "h.scores.tmp.step.""$STEP"".txt"
            mv "h.scores.tmp.step.""$STEP"".txt" "h.scores.step.""$STEP"".txt"

            # procede with LIGer: fix and overwrite confidence file, accept links and update assembly
            gawk -f $compute_confidences_script "h.scores.step.""$STEP"".txt" | sort -r -nk4 -nk5 > "h.double.sorted.confidence.step.""$STEP"".txt"
            gawk -f $accept_links_script "h.double.sorted.confidence.step.""$STEP"".txt" > "h.scaffolds.step.""$STEP"".txt"
            gawk -f $update_assembly_script "h.scaffolds.step.""$STEP"".txt" "h.scaffolds.original.notation.step.""$STEP"".txt" > "h.scaffolds.original.notation.tmp.step.""$STEP"".txt"
            mv "h.scaffolds.original.notation.tmp.step.""$STEP"".txt" "h.scaffolds.original.notation.step.""$STEP"".txt"

            # check if was helpful?
            if [[ -s "h.scaffolds.step.""$STEP"".txt" ]]; then
                STEP=$((1+$STEP))
                echo "...Starting iteration # $STEP"
                break
            fi

        done
    fi

done

# CONSOLIDATE FINAL OUTPUT
basenamefile="$(basename $contigPropFile .cprops)"
cp "h.scaffolds.original.notation.step.""$STEP"".txt" "$basenamefile"".asm"
for i in $(find . -maxdepth 1 -name 'h.dropouts.step.*.txt' | sort -t "." -nr -k 5); do
    cat $i >> "$basenamefile"".asm"
done


## CLEAN LEFTOVER HELPER FILES. TODO: delete some from inside the loop to save space.
find . -maxdepth 1 -name "h.*.txt" -delete

echo ":) DONE!"

##
##
##
##
##
##
##
