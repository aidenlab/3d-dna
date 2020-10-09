#### Description: Script to map 2D annotation files from input onto the assembly chromosome according to a particular assembly file, including non-diagonal cases
#### Usage: awk [ -v scale=<scale> ] [ -v ignore_id_orientation=1 ] -f lift-input-annotations-to-assembly-annotations.awk <assembly> <input-annotation-file>
#### Input: Assembly file and 2D annotation file (0-based, 12-column tab-separated assembly format).
#### Output: 2D annotation file STDOUT (0-based, 12 column tab-separated assembly format).
#### Parameters: Scale to stretch coordinates in order to fit into the assembly chromosome, can be both >1 (to squeeze the map) or <1 (to stretch the map out). Default scale = 1.
#### Notes: In principle can take in parameters for for fixed gap between individual contigs (cgap) and individual scaffolds (sgap), but not prompted. By default both are equal to 0.
#### WARNING: Works for diagonal annotations only at this point.
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu, version 20181007

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
(FILENAME==ARGV[1] && $0~/^>/){

	cname[substr($1,2)]=$2
	clength[$2]=$3
	next

}
#read in assembly
(FILENAME==ARGV[1]){

		for (i=1; i<=NF; i++)
		{	
			if($i<0){
				$i=-$i
				orientation[$i] = 16
			}else{
				orientation[$i] = 1
			}

			if (!( $i in clength ))
			{
				print ":( Assembly file does not match cprops file. Exiting!" > "/dev/stderr"
				print $0 > "/dev/stderr"
				exit 1
			}

			globalpos[$i] = pos

			pos += clength[$i]
			pos += cgap

		}
		pos += sgap
		next
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
	
	if (!($4 in cname))
	{
		next
	}
	else
	{
		$4=cname[$4]
	}
	
	if (!($1 in globalpos))
	{
		next
	}
	
	if (!($4 in globalpos))
	{
		next
	}

	# use unscaled coordinates, 1-based
	
	if($9==""){$9=$2}
	if($10==""){$10=$3}
	if($11==""){$11=$5}
	if($12==""){$12=$6}	

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

	}
	
	# use unscaled coordinates, 1-based
	$5=$11+1
	$6=$12

	if ( orientation[$4] == 1)
	{
		tmp5 = globalpos[$4] + $5 - 1
		tmp6 = globalpos[$4] + $6 - 1
	}
	else
	{
		tmp5=-$6 + globalpos[$4] + clength[$4]		
		tmp6=-$5 + globalpos[$4] + clength[$4]

	}

	## TODO: do I need this in non-diagonal case?
	if(!ignore_id_orientation && $4==$1 && orientation[$1] == 16){
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
	$9=tmp2-1; $10=tmp3; $11=tmp5-1; $12=tmp6;	# unscaled coordinates
	$2=int($9/scale); $3=int($10/scale);	$5=int($11/scale); $6=int($12/scale);	# scaled coordinates
	
	print
}

