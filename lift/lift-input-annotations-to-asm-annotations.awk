#### Description: Script to map 2D annotation files from input onto the assembly chromosome according to a particular asm and cprops file.
#### Usage: awk -v scale=<scale> -f remap-input-annotations-to-asm-annotations.awk <cprops> <asm> <input-annotation-file>
#### Input: Cprops file, asm file and 2D annotation file (0-based, 12-column tab-separated assembly format).
#### Output: 2D annotation file STDOUT (0-based, 12 column tab-separated assembly format).
#### Parameters: Scale to stretch coordinates in order to fit into the assembly chromosome, can be both >1 (to squeeze the map) or <1 (to stretch the map out). Default scale = 1.
#### Notes: In principle can take in parameters for for fixed gap between individual contigs (cgap) and individual scaffolds (sgap), but not prompted. By default both are equal to 0.
#### WARNING: Works for diagonal annotations only at this point.
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu, version 12/09/2016

BEGIN{
	OFS="\t"
	if (!scale)
	{
		scale=1
	}
#	cgap = 0
#	sgap = 0
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
#read in annotation file
{
	FS="\t"
	if (FNR==1)
	{
		print
		next
	}
	
	if (!($1 in cname))
	{
		next
	}
	else
	{
		$1=cname[$1]
	}
	
	if (!($1 in asm))
	{
		next
	}

	# use unscaled coordinates, 1-based
	$2=$9+1
	$3=$10

	if ( orientation[$1] == 1)
	{
		tmp2 = globalpos[$1] + $2 - 1
		tmp3 = globalpos[$1] + $3 - 1
	}
	else
	{
		tmp2=-$3 + globalpos[$1] + clength[$1]		
		tmp3=-$2 + globalpos[$1] + clength[$1]

		if ($8~/^-/)
		{
			sub("^\\-","+", $8)
			sub(" \\(\\-"," (+", $8)
		}
		else if ($8~/^+/)
		{
			sub("^\\+","-", $8)
			sub(" \\(\\+"," (-", $8)
		}
	}
	
	## TODO: Handle non-diagonal case!
	$1="assembly"; $4=$1;
	$9=tmp2-1; $10=tmp3; $11=$9; $12=$10;	# unscaled coordinates
	$2=int($9/scale); $3=int($10/scale);	$5=$2; $6=$3;	# scaled coordinates
	
	print
}

