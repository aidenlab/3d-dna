#### Description: Script to map mnd file onto the assembly chromosome according to a particular assembly file.
#### Usage: awk -v scale=<scale> -f lift-input-mnd-to-assembly-mnd.awk <assembly> <input-mnd>
#### Input: assembly file and mnd.
#### Output: Reduced assembly mnd stdout in standard 16 column format.
#### Parameters: scale to stretch or shrink coordinates (>1 to squeeze the map, <1 to stretch the map out). Default scale = 1. Stuff not fitting into the field of view (2.1 Gb) after scaling is ignored.
#### Notes: TODO: proper fragment number tracking for consistency.
#### Notes: In principle can take in parameters for for fixed gap between scaffolds, deprecated implementation.
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu on 20190102

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

#read in cprops part of assembly
( FILENAME == ARGV[1] ) && ( $0~/^>/ ){
	cname[substr($1,2)]=$2
	clength[$2]=$3
	next
}
#read in asm part of assembly
( FILENAME == ARGV[1] ){
	for (i=1; i<=NF; i++)
	{
		m=split($i, b, ",")
		shift=0
		for(j=1; j<=m; j++){
			
			if(!( abs(b[j]) in clength ))
			{
				print ":( Assembly file does not match cprops file. Missing "b[j]". Exiting!" > "/dev/stderr"
				exit 1
			}
			
			if( b[j]<0 ){
				b[j]=-b[j]
				orientation[b[j]] = -1
			}else{
				orientation[b[j]] = 1
			}
			globalpos[b[j]] = pos
			
			if (clength[b[j]]>shift){
				shift=clength[b[j]]
			}
		}

		
		pos += shift
		#pos += cgap

	}
	
	#pos += sgap
	next
}
{

	if (!($2 in cname) || !($6 in cname))
		next

	$2=cname[$2]
	$6=cname[$6]

	if (!($2 in orientation) || !($6 in orientation))
		next


	if (orientation[$2]==-1)
		$3 = (-$3 + globalpos[$2] + clength[$2])/scale
	else
		$3 = ($3 + globalpos[$2] - 1)/scale
	if (orientation[$6]==-1)
		$7 = (-$7 + globalpos[$6] + clength[$6])/scale
	else
		$7 = ($7 + globalpos[$6] - 1)/scale
		
	if ($3>CUTOFF || $7>CUTOFF){next}
	
	if ($2!=$6)
	{
		$4=0	## TODO: properly remap fragments?
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
}
1
