#### Description: Script to lift bed files files from input onto chromosome-length assembly described typically with _HiC.assembly file.
#### Usage: awk -v sandbox=<0|1> -v outlabel="label" -v editlabel="label" -f lift-input-bed-to-HiC-bed.awk <HiC.assembly> <input-bed-file>
#### Input: Chromosome-length assembly file and input bed file (can be 3-12 column).
#### Output: bed file (3-12 columns).
#### Parameters: Scale to stretch coordinates in order to fit into the assembly chromosome, can be both >1 (to squeeze the map) or <1 (to stretch the map out). Default scale = 1.
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu, version 10/15/2020

BEGIN{
	OFS="\t"
	pos = 1
	if (!outlabel){
		if(!sandbox)
			outlabel="assembly"
		else
			outlabel="HiC_scaffold"
	}
	if(!editlabel) editlabel=":::"
	if(!scale){scale=1}
}
# read in assembly file
FILENAME==ARGV[1]{
	if($0~/^>/){

		#read in cprops
		$1=substr($1,2)
		n=split($1, a, editlabel)
		clean_name=a[1]
		fragmentcounter[clean_name]+=1	# how many fragments are there corresponding to the original sequence
		fragmentid[clean_name" "fragmentcounter[clean_name]]=$2
		fragmentshift[clean_name" "fragmentcounter[clean_name]]=originallength[clean_name]
		originallength[clean_name]+=$3
		clength[$2]=$3
	} else {

		#read in assembly
		if(sandbox){
			chromcount++
			pos=1
		}

		for (i=1; i<=NF; i++) {
			if($i<0) {
				orientation[-$i]=-1
				$i=-$i
			} else {
				orientation[$i]=1
			}

			if (!( $i in clength )) {
				print ":( Assembly file does not match cprops file. Missing "$i". Exiting!" > "/dev/stderr" # should not happen
				exit 1
			}

			chrom[$i]=chromcount

			globalpos[$i] = pos

			pos += clength[$i]

		}
	}
	next
}
#read in bed file
{
	originalline=$0;

	FS="\t"
	if (!($1 in originallength)){
		# must be a comment line
		print ":| Warning: Interpreting the following as a comment line: " > "/dev/stderr"
		print originalline > "/dev/stderr"
		print; next
	}

	clean_name=$1

	if (fragmentcounter[clean_name]==1){
		# no splitting of original contig
		$1=fragmentid[clean_name" "1]
	} else {
		# contig was split, need to see where the interval falls

		i=1	# identify the fragment to which to assign the contig
		while (fragmentshift[clean_name" "(i+1)] < $2 && i < fragmentcounter[clean_name])
		{
			i++
		}

		j=1	# identify the fragment to which to assign the contig
		while (fragmentshift[clean_name" "(j+1)] < $3 && j < fragmentcounter[clean_name])
		{
			j++
		}

		if(i!=j)
		{
			print ":| Warning: interval was split. Skipping!" > "/dev/stderr"
			print originalline > "/dev/stderr"
			next
		}
		$1=fragmentid[clean_name" "i]
		$2-=fragmentshift[clean_name" "i]
		$3-=fragmentshift[clean_name" "i]
		if($6)
		{
			if(orientation[fragmentid[clean_name" "i]]<0)
			{
				if($6=="+")
				{
					$6="-"
				}else
					$6="+"
			}
		}
		if($7){$7-=fragmentshift[clean_name" "i]}
		if($8){$8-=fragmentshift[clean_name" "i]}
		if($12){k=split($12,a,","); $12=a[1]-fragmentshift[clean_name" "i]; for(j=2; j<=k; j++){if(a[j]!=""){$12=$12","(a[j]-fragmentshift[clean_name" "i])}}} #not tested
	}
	
	if (!($1 in chrom)){
		print ":| Warning: input fragment is not present in the output assembly:" > "/dev/stderr"
		print originalline > "/dev/stderr"
		next
	}

	# lift sequences to their final positions

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
	
	$2=int(tmp2/scale)
	$3=int(tmp3/scale)
	$1=outlabel"_"chrom[$1]
	
	print
}

