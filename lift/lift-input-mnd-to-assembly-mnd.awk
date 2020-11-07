#### Description: Script to map mnd file onto the assembly chromosome(s) according to a particular assembly file.
#### Usage: awk -v scale=<scale> -v sandbox=<0|1> -v label="label" -f lift-input-mnd-to-assembly-mnd.awk <assembly-file> <mnd-matching-assembly-file>
#### Input: Assembly file and mnd file matching assembly.
#### Output: Reduced assembly mnd stdout in standard 16 column format.
#### Parameters:
	#scale to stretch or shrink coordinates (>1 to squeeze the map, <1 to stretch the map out). Default scale = 1. Stuff not fitting into the field of view (2.1 Gb) after scaling is ignored. 				
	#global: true for lifting to assembly chromosome, false for lifting onto individual superscaffolds
	#label:	label for the chromosome(s). Default: HiC_scaffold.
#### Notes: TODO: proper fragment number tracking for consistency.
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu on 07/17/2016

function rcread(str,	a, i, n, outstr)
{
	n=split(str, a, "");
	for (i=n; i>0; i--)
	{
		if(a[i]=="A")
			outstr=outstr"T"
		else if (a[i]=="T")
			outstr=outstr"A"
		else if (a[i]=="C")
			outstr=outstr"G"
		else if (a[i]=="G")
			outstr=outstr"C"
		else
			outstr=outstr""a[i]
	}
	return outstr
}

function rccigar(str,	a, s, n, m, i, outstr)
{
	n = split(str, a, "[MIDNSHPX=]") # Only deal with M, I, D, H and S.
	m =	split(str, s, "[0-9]+")
		
	for(i=m; i>1; i--)
	{
		outstr=outstr""a[i-1]""s[i]
	}
	return outstr
}

BEGIN{
	if (!scale)
	{
		scale=1
	}
	CUTOFF=2100000000
	pos = 1
	if (!label && !sandbox)
	{
		label="assembly"
	}
}
function abs(value)
{
  return (value<0?-value:value)
}
#read in assembly
(FILENAME==ARGV[1]){
	
	# read in cprops part
	if($0~/^>/)
	{
		cname[substr($1,2)]=$2
		clength[$2]=$3
		next
	}

	# read in assembly part
	chromcount++
	
	if(sandbox) pos=1
	
	n=split($0,a)

	for (i=1; i<=n; i++)
	{
		if (!( abs(a[i]) in clength ))
		{
			print ":( Assembly file does not match cprops file. Exiting!" > "/dev/stderr"
			exit 1
		}

		chrom[abs(a[i])]=chromcount

		asm[abs(a[i])] = 1
		globalpos[abs(a[i])] = pos
		if (a[i]<0)
		{	
			orientation[abs(a[i])] = 16
		}
		else
		{
			orientation[abs(a[i])] = 1
		}	
		pos += clength[abs(a[i])]
	}
	next
}
#read in mnd file
{
	if (!($2 in cname) || !($6 in cname))
		next

	$2=cname[$2]
	$6=cname[$6]

	if (!($2 in asm) || !($6 in asm))
		next

	if (orientation[$2]==16)
	{
		$3 = (-$3 + globalpos[$2] + clength[$2])/scale
		if(verbose)
		{
			if ($1=="0"){$1="16"}
			else if ($1=="16"){$1="0"}
			$10=rccigar($10)
			$11=rcread($11)
		}
	}
	else
	{
		$3 = ($3 + globalpos[$2] - 1)/scale
	}
	if (orientation[$6]==16)
	{
		$7 = (-$7 + globalpos[$6] + clength[$6])/scale
		if(verbose)
		{
			if ($5=="0"){$5="16"}
			else if ($5=="16"){$5="0"}
			$13=rccigar($13)
			$14=rcread($14)
		}
	}
	else
	{
		$7 = ($7 + globalpos[$6] - 1)/scale
	}
	
	if (($3<=CUTOFF)&&($7<=CUTOFF)) # if scale is properly chosen this should always be true
	{
		if ($2!=$6)
		{
			$4=0	## TODO: maybe properly remap fragments
			$8=1
		}
		
		if(sandbox)
		{
			if(chrom[$2]<=chrom[$6])
			{
				$2=label""chrom[$2]
				$6=label""chrom[$6]
				$3 = int($3)
				$7 = int($7)
			}
			else
			{
				tmp=$2
				$2=label""chrom[$6]
				$6=label""chrom[tmp]
				tmp=$3
				$3 = int($7)
				$7 = int(tmp)			
			}
		}
		else
		{
			$2=label
			$6=label
			$3 = int($3+1)	#don't remember why +1?? because of scaling so that it's not 0?
			$7 = int($7+1)
		}

		#reduce size of output
		if (!verbose)
		{
			$10="-"
			$11="-"
			$13="-"
			$14="-"
			$15="-"
			$16="-"
		}
		print
	}
}
