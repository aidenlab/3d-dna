## Helper script to convert ps files to assembly files.
## USAGE: awk -f psf-to-assembly.awk <path_to_psf_file> 
## INPUT: psf file (">CHR POS REF ALT ID")
## OUTPUT: stdout in assembly format
## Written by O.D. (olga.dudchenko@bcm.edu)
## Version: 4/29/19

# read psf properties part (expected format ">CHR POS REF ALT ID")
$0~/^>/{
	print $1":"$2":"$3, 2*$NF-1, 1000
	print $1":"$2":"$4, 2*$NF, 1000
	next
}
# read psf set part
{
	ref=""
	alt=""
	for(i=1;i<=NF; i++)
	{
		if($i>0)
		{
			ref=ref" "(2*$i-1)
			alt=alt" "(2*$i)
		}
		else
		{
			ref=ref" "(-2*$i)
			alt=alt" "(-2*$i-1)
		}
	}
	print substr(ref,2)
	print substr(alt,2)
}