## Helper script to convert ps files to pseudo-vcf files.
## USAGE: awk -f psf-to-vcf.awk <path_to_psf_file> 
## INPUT: psf file
## OUTPUT: stdout in vcf format
## Written by O.D. (olga.dudchenko@bcm.edu)
## Version: 4/30/19
## NOTE: not real vcf in that there is no attached reference, variants are listed in arbitrary order in REF and ALT columns. Phase set id is position of first variant as per convention. All are assigned 100 quality score and PASS filter.

BEGIN{OFS="\t"}
# read psf properties part
$0~/^>/{
	chr[$5]=substr($1,2)
	pos[$5]=$2
	var1[$5]=$3
	var2[$5]=$4
	next
}
# read psf set part

## singletons
NF==1{print chr[$1], pos[$1], ".", var1[$1], var2[$1], 100, "PASS", ".", "GT", "0/1"; next}

## phased blocks
{
	if(pos[$1])
		ps=pos[$1]
	else	# should not happen, check
		ps=pos[-$1]
	for(i=1;i<=NF;i++)
	{
		if($i>0)
			print chr[$i], pos[$i], ".", var1[$i], var2[$i], 100, "PASS", ".", "GT:PS", "0|1:"ps
		else
			print chr[-$i], pos[-$i], ".", var1[-$i], var2[-$i], 100, "PASS", ".", "GT:PS", "1|0:"ps			
	}
}