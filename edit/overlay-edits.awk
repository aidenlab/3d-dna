#### Description: New version of the script that overlays mismatch_narrow track with input boundaries 
#### Usage: awk -v bin_size=<bin_size> -f overlay-edits.awk <input_annotation_file> <mismatch_narrow_bed>
#### Input: input_annotation file (asm coordinates), mismatch_narrow bed file (asm coordinates)
#### NOTE: the bed file is assumed to be merged (could merge here). bin_size is used to determine minimal acceptable edge overlap (assumed to be equal to bin_size of the mismatch bed file).
#### Output: stdout 2D annotation file with two types of labels: "debris" (for removal), "mismatch" (for invosertion analysis)
#### TODO: introduce different labels to distinguish between misassemblies and repeats
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu on 12/27/2016

function printout(a, b, flag){
	if (flag=="debris")
		color="255,255,0"
	else
		color="0,255,0"
	print "assembly", int(a/z), int(b/z), "assembly", int(a/z), int(b/z), color, flag, a, b, a, b 
}

BEGIN{
	FS="\t"
	OFS="\t"
	print "chr1", "sx1", "sx2", "chr2", "sy1", "sy2", "color", "id", "x1", "x2", "y1", "y2"
}

# read in annotations
FILENAME==ARGV[1]{
	if(FNR==1){next}
	
	if(FNR==2){z=int($10/$3)}
	if(($2!=0 && int($9/$2)!=z) || z==0)	# should not happen
	{
		print ":( Cannot interpret scaling in the annotation file. Exiting! " > "/dev/stderr"
		exit
	}
	
	start_scaf[FNR-2]=$9 # use 'actual' coordinates
	end_scaf[FNR-2]=$10
	total_scaffolds=FNR-1
	next
}
FNR==1{asort(start_scaf); asort(end_scaf)}
# read in the bed file
{
		# use 'actual' coordinates
		$2=z*$2
		$3=z*$3
		
		# find scaffold in which lies the start of the depleted region
		while (end_scaf[i] <= $2)
		{
			i++
		}
		if ($3<=end_scaf[i])	# intra-input scaffold depletion, flag for removal
		{
			printout($2, $3, "debris")
#			print $2, $3, start_scaf[i], end_scaf[i], "intra-scaffold"
		}
		else	# inter-input scaffold depletion
		{		
			# check if overlap with current scaffold is long
			if (end_scaf[i]-$2 > z*bin_size)	# might be causative, flag for removal
			{
				## NOTE: could have checked if leftover is too small and remove without breaking
				printout($2, end_scaf[i], "debris")
			}
			else
			{
				printout($2, end_scaf[i], "mismatch")
			}
			# scaffolds that fall into the annotated region
			i+=1
			while (end_scaf[i] < $3 && i < total_scaffolds)
			{
				printout(start_scaf[i], end_scaf[i], "debris")
				i++
			}

			if ($3-start_scaf[i] > z*bin_size)
			{
				## NOTE: could have checked if leftover is too small and remove without breaking
				printout(start_scaf[i], $3, "debris")
			}
			else
			{
				printout(start_scaf[i], $3, "mismatch")
			}
		}
}
