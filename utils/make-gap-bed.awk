#### Description: Awk script to annotate gaps in fasta.
#### USAGE: awk -f make-gap-bed.awk <path-to-fasta>
#### Options: can take in gap_size_threshold. Default is 0 (all Ns or ns are annotated).
#### Input: fasta file
#### Output: bed-formatted stdout
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu.
BEGIN{
#	gap_size_threshold=500
}
{
	if ($0~/>/)	# new scaffold
	{
		if ((start)&&(counter >= gap_size_threshold))	# in big gap
			print scaf_id, start-1, start+counter-1 # bed start is 0-based and bed end is 1-based

		scaf_id=substr($1, 2)
		pos=0	# position in contig
		start=0	# gap start
		counter=0	# gap length
		next
	}
	
	if ($0!~/N/&&$0!~/n/)
	{
		if ((start)&&(counter>=gap_size_threshold))
			print scaf_id, start-1, start+counter-1
		start=0		
		counter=0
		pos+=length($0)
	}
	else
	{
		n=split($0, a, "")
		for (i=1; i<=n; i++)
		{
			if (a[i]=="N" || a[i]=="n")
			{
				if (!start)
				{
					start=pos+1
				}
				counter+=1
			}
			else
			{
				if ((start)&&(counter>=gap_size_threshold))
					print scaf_id, start-1, start+counter-1
				start=0			
				counter=0
			}
			pos+=1
		}
	}	
}
END{
	if ((start)&&(counter>=gap_size_threshold))
		print scaf_id, start-1, start+counter-1
}
