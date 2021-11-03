## single-thread script
## USAGE: awk [ -v mapq=1 ] -f extract-SNP-reads-from-sam-file.awk <path-to-psf-file> <path-to-sam-file>
## INPUT: psf file in ">chr pos var1 var2 id" format. IMPORTANT: psf file is assumed to be sorted by chr and position.
## INPUT: sam file or pipe
## OUTPUT: 
## OPTIONS: mapq [default: 1]. Parameter for minimal read mapping quality threshold.
## OD, original version: 210817.
## NOTE: Chromosome names should not contain ":"!

BEGIN{
	# defaults
	if(length(mapq) == 0){mapq=1}
	# other params are 0 by default
}

function shift(ARRAY,	i, n)
{
	for (i in ARRAY)
		n++
	for (i=1; i<=n-1; i++)
	{
		ARRAY[i] = ARRAY[i+1]
	}
	delete ARRAY[i]
}

# read psf file properties section
FILENAME==ARGV[1]{
	if($0~/^>/)
	{
		$1=substr($1,2)":"$2	# for now, rewrite later
		$2=$3
		$3=$4
	
		split($1, a, ":")
		counter[a[1]]++
		snp[a[1]" "counter[a[1]]] = $1
		position[$1] = a[2]
		if(length($2)>=length($3))
		{
			var1[$1] = $2
			var2[$1] = $3
		}
		else
		{
			var1[$1] = $3
			var2[$1] = $2	
		}
	}
	next
}

# read sam file

{
	id_field=$1
	flag_field=$2
	chr_field=$3
	pos_field=$4
	mapq_field=$5
	cigar_field=$6
	sequence_field=$10

	
	if((!chr_field in counter)||mapq_field<mapq){next}
	
	# Ugly workaround for old gawk, split does not spit out separators, have to split separately and shift
	cigarlength = split(cigar_field, cigar, "[MIDNSHPX=]") # Only deal with M, I, D, H and S.
	split(cigar_field, seps, "[0-9]+")		
	shift(seps)	
	cigarlength-=1

	# ignore pairs with weird cigars && count total reference match length
	matchlength=0;
	for(i in seps)
	{
		if(seps[i]!="M"&&seps[i]!="S"&&seps[i]!="I"&&seps[i]!="D"&&seps[i]!="H")
		{
			next	# this should not happen.
		}
		if (seps[i] == "M" || seps[i] == "D")
			matchlength += cigar[i]
	}

	# get leftmost and rightmost match positions on the reference from mnd file, ignore weird flags
	leftmatch=pos_field
	rightmatch=leftmatch+matchlength-1

	# half-search for at least one overlapping snp
	upstream=1
	downstream=counter[chr_field]
	
	if(leftmatch > position[snp[chr_field" "downstream]] || rightmatch < position[snp[chr_field" "upstream]])
	{
		next
	}

	while (1)
	{
		tmp = int ((upstream+downstream)/2)
		if (rightmatch < position[snp[chr_field" "tmp]])
			downstream = tmp
		else if (leftmatch > position[snp[chr_field" "tmp]])
			upstream = tmp
		else
			break

		no_overlaps=0;
		# no overlapping snps
		if (upstream+1 >= downstream)
		{
			no_overlaps=1;
			break;
		}
	}

	if (no_overlaps)
	{
		next
	}

	# track back in case this is not the first overlapping snp
	while (leftmatch <= position[snp[chr_field" "(tmp-1)]] && tmp > 1) {
		tmp--
	}

	# parse sequence to see which variant(s) show up in sequence	
	snplist=""

	# string position tracker
	s=1 
	for (k=1; k<=cigarlength; k++)
	{
		if (seps[k] == "S")
		{
			s += cigar[k]
		}
		else if (seps[k] == "H")
		{
			continue
		}
		else if (seps[k] == "D")
		{
			leftmatch += cigar[k]
		}
		else if (seps[k] == "I")
		{
			s += cigar[k]
		}
		else
		{
			## the whole length match is not relevant for SNP comparison. Should I keep for the same of other vars or give up?
			if (position[snp[chr_field" "tmp]]<=leftmatch+cigar[k]-1)
			{
				if(substr(sequence_field, s + (position[snp[chr_field" "tmp]] - leftmatch), length(var1[snp[chr_field" "tmp]]))==var1[snp[chr_field" "tmp]])
				{
					snplist=snplist" "snp[chr_field" "tmp]":"var1[snp[chr_field" "tmp]]
				}
				else if (substr(sequence_field, s + (position[snp[chr_field" "tmp]] - leftmatch), length(var2[snp[chr_field" "tmp]]))==var2[snp[chr_field" "tmp]])
				{
					snplist=snplist" "snp[chr_field" "tmp]":"var2[snp[chr_field" "tmp]]
				}
			
				## look for more overlapping snps
				if (position[snp[chr_field" "(tmp+1)]]<=rightmatch && tmp<counter[chr_field])
				{
					tmp++
					k--
					continue
				} else {
					break					
				}
			}
			leftmatch+=cigar[k]
			s+=cigar[k]
		}
	}

	## end of read handling, report part (may change)

	n=split(substr(snplist,2),a," ")
	
	if(n==0){next}

	## report once with SNP randomly chosen from the overlap set. Might change this to something compatible with combinatorial reporting and/or amplification of signal from signle-alignment reads [have not been using in practice lately].   
	n=int(1+rand()*n)
	$3=a[n]
	print $0

}
