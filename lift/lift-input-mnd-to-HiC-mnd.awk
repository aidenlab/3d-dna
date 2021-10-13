#### Description: Script to edit and lift mnd files from input onto assembly based on the assembly file.
#### Usage: awk -v sandbox=<0|1> -v longformat=<0|1> -v scale=<?> -v mapq=<mapq> -f lift-input-mnd-to-HiC-mnd.awk <assembly-file> <input-mnd-file>
#### Input: assembly file and input mnd file (long format).

#### Output: mnd file (16 columns by default, 8 if short is invoked).
#### Parameters (prompted): sandbox. Default sandbox = 0 (assembly coordinates in the output).
#### Parameters (prompted): longformat=<0|1>. If true long format for output mnd is used. Default: 0, i.e. short format, with 8 columns, isused for output.
#### Parameters (prompted): Scale to stretch coordinates in order to fit into the assembly chromosome, can be both >1 (to squeeze the map) or <1 (to stretch the map out). Default scale = 1.
#### Unprompted: -v chromlabel="label" [default: "assembly" for sandbox=0 and "HiC_scaffold_" for sandbox=1] -v editlabel="label" [default: ":::"] -v verbose=<0|1> [default: 0].
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu, version 10/15/2020
#### WARINING: shold formally set separator for the mnd file, but 99.9999% of the time should be fine.
#### TODO: potentially could add proper rc for read sequences and cigars for long format. Some of it is implemented in lift-input-mnd-o-assembly-mnd.awk

BEGIN{
	pos = 1

	if(sandbox)
	{
		if (!chromlabel)
			chromlabel="HiC_scaffold_"
		else if (chromlabel=="null")
			chromlabel=""
	}else
	{
		if (chromlabel!=""&&chromlabel!="assembly")
			print ":| Warning: chromlabel setting is not compatible when not in sandbox output mode. Switching to default \"assembly\" label." > "/dev/stderr"
			
		chromlabel="assembly"
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
#read in mnd file
($9<mapq||$12<mapq){next}
{
	outputline=""

	if (!($2 in originallength) || !($6 in originallength)){
		if (verbose)
		{
			print ":| Warning: unexpected sequence name found in the mnd file. Skipping line:" > "/dev/stderr"
			print $0 > "/dev/stderr"
		}
		next
	}

	if($4!=$8){$4=0; $8=1} #TODO: maybe do proper fragment remapping

	for (half=1; half<=2; half++)
	{
		if (half==1)
		{
			orientation_field=$1
			chr_field=$2
			pos_field=$3
			fragment_field=$4
			mapq_field=$9
			cigar_field=$10
			sequence_field=$11
			readid_field=$15
		}
		else
		{
			orientation_field=$5
			chr_field=$6
			pos_field=$7
			fragment_field=$8
			mapq_field=$12
			cigar_field=$13
			sequence_field=$14
			readid_field=$16
		}

		if (pos_field>originallength[chr_field]){pos_field=originallength[chr_field]} #cigar shift correction

		clean_name=chr_field

		if (fragmentcounter[clean_name]==1)
		{
			# no splitting of original contig
			chr_field=fragmentid[clean_name" "1]
		} else
		{
			# contig was split, need to see where the interval falls

			i=1	# identify the fragment to which to assign the contig
			while (fragmentshift[clean_name" "(i+1)] < pos_field && i < fragmentcounter[clean_name])
			{
				i++
			}
			# perhpas revisit: may need more careful mapping of the position and skipping of split reads
			# j=1	# identify the fragment to which to assign the contig
			# while (fragmentshift[clean_name" "(j+1)] < $3 && j < fragmentcounter[clean_name])
			# {
			# 	j++
			# }

			# if(i!=j)
			# {
			# 	print ":| Warning: interval was split. Skipping!" > "/dev/stderr"
			# 	print originalline > "/dev/stderr"
			# 	next
			# }
			chr_field=fragmentid[clean_name" "i]
			pos_field-=fragmentshift[clean_name" "i]
			if(orientation[fragmentid[clean_name" "i]]<0)
			{
				if(orientation_field==0)
				{
					orientation_field=16
				}else
					orientation_field=0
			}
		}
	
		if (!(chr_field in chrom))
		{
			if(verbose)
			{
				print ":| Warning: input fragment is not found in the output assembly:" > "/dev/stderr"
				print $0 > "/dev/stderr"
			}
			next
		}

		# lift sequences to their final positions

		if ( orientation[chr_field] == 1)
		{
			pos_field = globalpos[chr_field] + pos_field - 1
		}
		else
		{
			pos_field=-pos_field + globalpos[chr_field] + clength[chr_field]		
		}
		pos_field=int(pos_field/scale)

		if(!sandbox)
			chr_field=chromlabel
		else
			chr_field=chromlabel""chrom[chr_field]

		outputline=outputline" "orientation_field" "chr_field" "pos_field" "fragment_field
	}

	if (longformat)
		outputline=outputline" "$9" "$10" "$11" "$12" "$13" "$14" "$15" "$16

	print substr(outputline, 2)
}

