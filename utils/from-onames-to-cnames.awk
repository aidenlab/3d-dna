#### Description: Helper script to switch from external contig/scaffold names to internal names.
#### USAGE: awk -f from-onames-to-cnames.awk <path-to-cprops> <path-to-asm>
#### Input: cprops and asm files
#### Output: asm-formatted stdout
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu.
FILENAME==ARGV[1]{
	cname[$1]=$2
	next
}
{
	str=""
	for(i=1; i<=NF; i++)
	{
		pref=""
		if ($i~/^-/)
		{
			pref="-"
			$i=substr($i,2)
		}
		str=str" "pref""cname[$i]
	}
	print substr(str,2)
}
