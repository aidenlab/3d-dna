#### Description: Script to map 2D annotation files from a given assembly to the original input scaffolds/contigs according to a particular asm and cprops file.
#### Usage: awk -v scale=<scale> -f remap-asm-annotations-to-input-annotations.awk <cprops> <asm> <asm-annotation-file>
#### Input: Cprops file, asm file and 2D annotation file. Note that annotations are expected to be in the bed-like 0-based coordinate format (0-1000; 1000-2000; 2000-3000), 12-column assembly tab separated format.
#### WARNING: Works for diagonal annotations only at this point.
#### WARNING: Annotations are expected not to cross input scaffold/contig boundaries. TODO: add splitting by boundary into this script rather than handle it outside.
#### Output: 2D annotation file STDOUT in assembly 12-column tab separated format.
#### Parameters: parameter to scale coordinates in final output, can be both >1 (to squeeze the map) or <1 (to stretch the map out, default scale = 1). Never needed, kept just in case for cases when individual scaffolds are > 2.1 Gb.
#### Notes: In principle can take in parameters for fixed gap between individual contigs (cgap) and individual scaffolds (sgap), but not prompted. By default both are equal to 0.
#### TODO: Maybe change globalpos from being defined as a function of cprops entry to function of counter, safer.. 
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu, version date 12/09/2016

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
#		cname[$1] = $2
		oname[$2] = $1
		clength[$2] = $3
		next
	}
}
#read in asm
{
	if (FILENAME==ARGV[2])
	{
		n=split($0,a)

		for (i=1; i<=n; i++)
		{
			if (!( abs(a[i]) in clength ))
			{
				print ":( Assembly file does not match cprops file. Exiting!" > "/dev/stderr"
				print $0 > "/dev/stderr"
				exit 1
			}
			counter += 1
			asm[counter] = abs(a[i])
			if ( abs(a[i]) in globalpos )
			{
				print ":( Assembly file contains multiple entries of the same contig. Exiting!" > "/dev/stderr"
				exit 1
			}
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
	
	# lift unscaled coordinates, switch to 1-based coordinate system
	$2=$9+1
	$3=$10
	
	# figure out which contig. Maybe TODO: optimize by sorting the original annotation file

	i = 1
	while ( globalpos[asm[i]] + clength[asm[i]] -1 < $2 && i <= counter)
	{
		i++
	}
	
	if ( i > counter)
	{
		print ":( Annotations are our of range covered by suggested contigs. Exiting!" > "/dev/stderr"
		exit 1
	}

	# TODO: handle splitting of annotations here.	
	if ( $3 > globalpos[asm[i]] + clength[asm[i]] )
	{
		print ":( Annotations span several contigs. Exiting!" > "/dev/stderr"
		print $0 > "/dev/stderr"
		
		print i, asm[i], clength[asm[i]], globalpos[asm[i]]
		
		exit 1
	}
	
	$1 = asm[i]

	if ( orientation[$1] == 1)
	{
		tmp2 = $2 - globalpos[$1] + 1
		tmp3 = $3 - globalpos[$1] + 1
		
	}
	else
	{
		tmp2 = -$3 + globalpos[$1] + clength[$1]
		tmp3 = -$2 + globalpos[$1] + clength[$1]

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
	
	#Diagonal case only. TODO: handle off-diagonal case!
	
	$1=oname[$1]; $4=$1;
	$9=tmp2-1; $10=tmp3; $11=$9; $12=$10;		# unscaled coordinates
	$2=int($9/scale); $3=int($10/scale);	$5=$2; $6=$3;	# scaled coordinates
	
	print
}
