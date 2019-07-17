#### Description: Script to map mnd file onto the assembly chromosome according to a particular asm file.
#### Usage: awk -v scale=<scale> -f lift-contig-mnd-to-asm-mnd.awk <cprops> <asm> <contig-mnd>
#### Input: Cprops file, asm file and mnd.
#### Output: Reduced assembly mnd stdout in standard 16 column format.
#### Parameters: scale to stretch or shrink coordinates (>1 to squeeze the map, <1 to stretch the map out). Default scale = 1. Stuff not fitting into the field of view (2.1 Gb) after scaling is ignored.
#### Notes: TODO: proper fragment number tracking for consistency.
#### Notes: In principle can take in parameters for for fixed gap between scaffolds, deprecated implementation.
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu on 07/17/2016

BEGIN{
	if (!scale)
	{
		scale=1
	}
	CUTOFF=2100000000
	pos = 1
}
function abs(value)
{
  return (value<0?-value:value)
}

#read in cprops
{
	if (FILENAME==ARGV[1])
	{
		cname[$1]=$2
		clength[$2]=$3
		next
	}
}
#read in assembly
{
	if (FILENAME==ARGV[2])
	{
		n=split($0,a)

		for (i=1; i<=n; i++)
		{
			if (!( abs(a[i]) in clength ))
			{
				print ":( Assembly file does not match cprops file. Exiting!" > "/dev/stderr"
				exit 1
			}

			asm[abs(a[i])] = 1
			globalpos[abs(a[i])] = pos
			if (abs(a[i])!=a[i])
			{	
				orientation[abs(a[i])] = 16
			}
			else
			{
				orientation[abs(a[i])] = 1
			}	
			pos += clength[abs(a[i])]
			pos += cgap

		}
		pos += sgap
		next
	}
}
{
	if (!($2 in cname) || !($6 in cname))
		next

	$2=cname[$2]
	$6=cname[$6]

	if (!($2 in asm) || !($6 in asm))
		next

	if (orientation[$2]==16)
		$3 = (-$3 + globalpos[$2] + clength[$2])/scale
	else
		$3 = ($3 + globalpos[$2] - 1)/scale
	if (orientation[$6]==16)
		$7 = (-$7 + globalpos[$6] + clength[$6])/scale
	else
		$7 = ($7 + globalpos[$6] - 1)/scale
	
	if (($3<=CUTOFF)&&($7<=CUTOFF)) # if scale is properly chosen this should always be true
	{
		if ($2!=$6)
		{
			$4=0	## TODO: properly remap fragments
			$8=1
		}
		$2="assembly"
		$6="assembly"

		$3 = int($3+1)
		$7 = int($7+1)

		#reduce size of output
		$11="R1"
		$14="R2"
		$15="Id1"
		$16="Id2"
		print $0
	}
}
