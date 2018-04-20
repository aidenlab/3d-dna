#!/bin/bash

## Wrapper script to merge contigs inside clusters as indicated by the annotated asm file.
## Requires LASTZ
## Written by: OD

USAGE="
*****************************************************
USAGE: merge-tiled-asm.sh -a <tiled-annotations> <path-to-cprops> <path-to-annotated-asm> <path-to-fasta-split-folder>
*****************************************************
"
# DEPENDENCIES: handle better!
pipeline=`cd "$( dirname $0)" && cd .. && pwd`

# DEFAULTS:
merger_lastz_options=\"--gfextend\ --gapped\ --chain=200,200\"

# coloring for the annotation file, probably overkill
cc_intra="255,255,255"
cc_break="255,0,0"
cc_overlap="255,255,0"

## HANDLE OPTIONS
while getopts "a:o:h" opt; do
case $opt in
	h) echo "$USAGE"
		exit 0
	;;
	a) if [ -f $OPTARG ]; then
			echo ":) -a flag was triggered, using $OPTARG annotations to build a qc annotation track" >&1
			tiled_annotations=`readlink -f $OPTARG`
		else
			echo ":( Could not find the file $OPTARG. Will skip building a qc annotation track" >&2
		fi
	;;
	o) re='^\"--.+\"$'
		if [[ $OPTARG =~ $re ]]; then
			echo ":) -o flag was triggered, assuming ${OPTARG} as a list of options to pass on to LASTZ to tune alignment." >&1
			merger_lastz_options=${OPTARG}
		else
			echo ":( Wrong syntax for LASTZ option string. Using the default value ${merger_lastz_options}." >&2
		fi
	;;
	*) echo "$USAGE" >&2
		exit 1
	;;
esac
done
shift $(( OPTIND-1 ))

# remove comments from LASTZ options
merger_lastz_options="${merger_lastz_options%\"}"
merger_lastz_options="${merger_lastz_options#\"}"

## HANDLE ARGUMENTS
[ -z $1 ] || [ -z $2 ] || [ -z $3 ] && echo >&2 ":( Some input seems missing." && echo >&2 "$USAGE" && exit 1

cprops=`readlink -f $1`
asm=`readlink -f $2`
splitFas=`readlink -f $3`

filename=`basename $cprops .cprops`

## PREP WORKING DIR
chrname=${filename}
outDir="${chrname}"
[ -d $outDir ] && rm -r ${outDir}
mkdir $outDir && cd $outDir

main_fa_file="merged_${chrname}.fa"
tmp_merge_file="tmp_merged_${chrname}.fa"
tmp_fa_file="tmp_${chrname}.fa"
merged_asm=${chrname}"_merged.asm"
touch "${main_fa_file}"

echo "chr1	x1	x2	chr2	y1	y2	color	id	X1	X2	Y1	Y2" > shifts_2D_input.txt

merge="false"
contig_counter=0
while read -r line
do
#echo $line
	first=$((contig_counter+1))
	for contig in $line
	do
#echo $contig
		anno_contig=$contig
		if [[ $contig = \{* ]]; then
			contig=$(echo $contig | cut -d'{' -f 2)
		fi
		if [[ $contig = *\} ]]; then
			contig=$(echo $contig | rev | cut -d'}' -f 2 | rev)
		fi

#echo $contig
		## might use more annotations
		if [[ $contig == -* ]]; then
			reverse=1
			contig=$(echo $contig | cut -d'-' -f 2)
		else
			reverse=0
		fi
		
		
		if [ ${contig_counter} -eq 0 ]; then	## handle first
			
			if [ $reverse -eq 0 ]; then
				cat ${splitFas}/${contig}.fa > "${tmp_merge_file}"
			else
				awk -f ${pipeline}/utils/reverse-fasta.awk ${splitFas}/${contig}.fa > "${tmp_merge_file}"
			fi
			globalpos=0
			len=$(awk -f ${pipeline}/utils/generate-cprops-file.awk ${splitFas}/${contig}.fa | awk '{print $NF}')
			let contig_counter=contig_counter+1
			
			## TODO: redo this via shifts, scaling also needs to be incorporated
			echo "assembly	$(( globalpos ))	$(( globalpos+len ))	assembly	$(( globalpos ))	$(( globalpos+len ))	${cc_break}	${contig}	$(( globalpos ))	$(( globalpos+len ))	$(( globalpos ))	$(( globalpos+len ))"	 >> "shifts_2D_input.txt"
			globalpos=$(( globalpos+len ))	
			last_match_pos=1
			last_tiled_pos=1	
		else
			if [ $reverse -eq 0 ]; then
				contig_file=${splitFas}/${contig}.fa
			else
				awk -f ${pipeline}/utils/reverse-fasta.awk ${splitFas}/${contig}.fa > "RC_contig.fa"
				contig_file="RC_contig.fa"
			fi
#head ${contig_file}			
			if [ $merge = "true" ]; then	# attempt to tile next if in cluster				
				align=$(lastz "${tmp_merge_file}" ${contig_file} ${merger_lastz_options} | awk -v last_tiled_pos=${last_tiled_pos} -f ${pipeline}/merge/extract-highest-oriented-tiled-stanza.awk)				
			else
				align=""
			fi
#echo $align >> alignments.txt
			read ts te prevlen qs qe len score <<< "$align"
			if [ -z "$ts" ]; then	# not merged or no overlap found
				echo "...No reliable overlap found. Breaking cluster."
				awk -f ${pipeline}/utils/wrap-fasta-sequence.awk "${tmp_merge_file}" >> "${main_fa_file}"
				if [ $reverse -eq 0 ]; then
					cat ${splitFas}/${contig}.fa > "${tmp_merge_file}"
				else
					awk -f ${pipeline}/utils/reverse-fasta.awk ${splitFas}/${contig}.fa > "${tmp_merge_file}"
				fi
				last_match_pos=1
				last_tiled_pos=1
				### add annotations
				len=$(awk -f ${pipeline}/utils/generate-cprops-file.awk ${splitFas}/${contig}.fa | awk '{print $NF}')
				let contig_counter=contig_counter+1
				echo "assembly	$(( globalpos ))	$(( globalpos+len ))	assembly	$(( globalpos ))	$(( globalpos+len ))	${cc_break}	${contig}	$(( globalpos ))	$(( globalpos+len ))	$(( globalpos ))	$(( globalpos+len ))" >> "shifts_2D_input.txt"
#				echo ${contig}"	" >> ${chrname}"_new_annotations.txt"
				globalpos=$(( globalpos+len ))
			else	# merging
				echo "...Overlap found!"
				last_match_pos=${ts}
				last_tiled_pos=$((ts-qs))
				
#echo $align
#echo  $(( len-prevlen+te-qe ))
				
				if [ $(( len-prevlen+te-qe )) -le 0 ]; then
					echo "....new scaffold overlaps fully with previous ones. Skipping merge."
					# new contig fully inside
					echo "assembly	$(( globalpos-prevlen+ts-qs ))	$(( globalpos+len-prevlen+te-qe ))	assembly	$(( globalpos-prevlen+ts-qs ))	$(( globalpos+len-prevlen+te-qe ))	${cc_intra}	${contig}	$(( globalpos-prevlen+ts-qs ))	$(( globalpos+len-prevlen+te-qe ))	$(( globalpos-prevlen+ts-qs ))	$(( globalpos+len-prevlen+te-qe ))" >> "shifts_2D_input.txt"
				else
					### second contig needs to be incorporated
					echo "....merging scaffolds"

					echo "assembly	$(( globalpos-prevlen+ts-qs ))	$(( globalpos+len-prevlen+te-qe ))	assembly	$(( globalpos-prevlen+ts-qs ))	$(( globalpos+len-prevlen+te-qe ))	${cc_overlap}	${contig}	$(( globalpos-prevlen+ts-qs ))	$(( globalpos+len-prevlen+te-qe ))	$(( globalpos-prevlen+ts-qs ))	$(( globalpos+len-prevlen+te-qe ))" >> "shifts_2D_input.txt"
					awk -v test=1 -v pos=${te} -f ${pipeline}/merge/grab-fix-sequence.awk ${tmp_merge_file} > ${tmp_fa_file}
					awk -v test=0 -v pos=$((qe+1)) -f ${pipeline}/merge/grab-fix-sequence.awk ${contig_file} | awk 'NR>1' >> ${tmp_fa_file}
					awk -f ${pipeline}/utils/wrap-fasta-sequence.awk ${tmp_fa_file} > ${tmp_merge_file}
					globalpos=$(( globalpos+len-prevlen+te-qe ))
				fi

			fi	


		fi
		
		if [[ $anno_contig = \{* ]]; then
			merge="true"
		fi
		if [[ $anno_contig = *\} ]]; then
			merge="false"
		fi
	done
	seq ${first} ${contig_counter} | xargs >> ${merged_asm}
done < $asm

awk -f ${pipeline}/utils/wrap-fasta-sequence.awk ${tmp_merge_file} >> ${main_fa_file}

[ -z ${tiled_annotations} ] || (awk -f ${pipeline}/merge/lift-merged-annotations-to-unmerged-map-shifts.awk ${tiled_annotations} shifts_2D_input.txt > shift_qc_track_2D_asm.txt && awk -f ${pipeline}/merge/lift-merged-annotations-to-unmerged-map-overlaps.awk ${tiled_annotations} shifts_2D_input.txt > overlap_qc_track_2D_asm.txt)

[ -f ${tmp_merge_file} ] && rm ${tmp_merge_file}
[ -f ${tmp_fa_file} ] && rm ${tmp_fa_file}
