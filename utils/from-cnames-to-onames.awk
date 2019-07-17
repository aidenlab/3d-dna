#### Description: Helper script to switch from external contig/scaffold names to internal names.
#### USAGE: awk -f from-cnames-to-onames.awk <path-to-cprops> <path-to-asm>
#### Input: cprops and asm files
#### Output: asm-formatted stdout
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu.
FILENAME==ARGV[1]{
	oname[$2]=$1
	next
}
{
	str=""
	for(i=1; i<=NF; i++)
	{
		if ($i~/^-/)
		{
			pref="-"
			$i=substr($i,2)
		}
		else
			pref=""
		if(oname[$i]!="")
			str=str" "pref""oname[$i]
	}
	print substr(str,2)
}
