## single-thread script
## USAGE: awk [ -v mapq=1 ] [ -v keep_native_coordinates=0/1 ] [ -v report_combinatorial_edges=0/1 ] [ -v amplify_neighboring_snps=0/1 ] [ -v keep_dangling_edges=0/1 ] -f extract-SNP-reads-from-mnd-file.awk <path-to-psf-file> <path-to-mnd-file>
## INPUT: psf file in ">chr pos var1 var2 id" format. IMPORTANT: psf file is assumed to be sorted by chr and position.
## INPUT: mnd file (16-format)
## OUTPUT: mnd stdout with chr calls substituted for snp ids in "chr:pos:var" format
## OPTIONS: mapq [default: 1]. Parameter for minimal read mapping quality threshold.
## OPTIONS: keep_native_coordinates [default: 0]. If false the coordinates are changed to 1 for downstread SNP-based contact mapping. If true coordinates are kept for contact mapping in 'native' coordinates.
## OPTIONS: report_combinatorial_edges [default: 0]. If false each read will correspond to at most one SNP (randomly chosen), even if the read overlaps with more than one SNP.
## OPTIONS: amplify_neighboring_snps [default: 0]. If true two edges will be added for SNPs that live on the same read/alignment (as opposed to reads observed on two different reads/alignments).
## OPTIONS: keep_dangling_edges [default: 0]. If true read pairs in which only one of the two reads overlap a SNP will be printed (with original chromosome name in a corresponding field).
## OD, original version: 2014; modified: 190417; latest: 210125.
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

# read mnd file

{	
	# handle read halves
	for (half=1; half<=2; half++)
	{
		if (half==1)
		{
			orientation_field=$1
			chr_field=$2
			pos_field=$3
			mapq_field=$9
			cigar_field=$10
			sequence_field=$11
		}
		else
		{
			orientation_field=$5
			chr_field=$6
			pos_field=$7
			mapq_field=$12
			cigar_field=$13
			sequence_field=$14
		}
	
	
		snplist[half]=""

		if((!chr_field in counter)||mapq_field<mapq){continue}
	
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

		if (orientation_field==0)
		{
			leftmatch=pos_field
			if (seps[1]=="S"||seps[1]=="H") leftmatch+=cigar[1]
			rightmatch=leftmatch+matchlength-1
		}
		else if (orientation_field==16)
		{
			rightmatch=pos_field
			if(seps[cigarlength]=="S"||seps[cigarlength]=="H") rightmatch-=cigar[cigarlength]
			leftmatch=rightmatch-matchlength+1
		}else
		{
			next	# this should not happen
		}
	# 	half-search for at least one overlapping snp

		upstream=1
		downstream=counter[chr_field]
	
		if(leftmatch > position[snp[chr_field" "downstream]] || rightmatch < position[snp[chr_field" "upstream]])
		{
			continue
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
			continue
		}
	
	# track back in case this is not the first overlapping snp

		while (leftmatch <= position[snp[chr_field" "(tmp-1)]] && tmp > 1) {
			tmp--
		}

	# 	print $0
	#	print snp[chr_field" "tmp], leftmatch, rightmatch

	# parse sequence to see which variant
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

				if (position[snp[chr_field" "tmp]]<=leftmatch+cigar[k]-1)
				{
					if(substr(sequence_field, s + (position[snp[chr_field" "tmp]] - leftmatch), length(var1[snp[chr_field" "tmp]]))==var1[snp[chr_field" "tmp]])
					{
						snplist[half]=snplist[half]" "snp[chr_field" "tmp]":"var1[snp[chr_field" "tmp]]
					}
					else if (substr(sequence_field, s + (position[snp[chr_field" "tmp]] - leftmatch), length(var2[snp[chr_field" "tmp]]))==var2[snp[chr_field" "tmp]])
					{
						snplist[half]=snplist[half]" "snp[chr_field" "tmp]":"var2[snp[chr_field" "tmp]]
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
	}

	n=split(substr(snplist[1],2),a," ")
	m=split(substr(snplist[2],2),b," ")
	
	if(n==0&&m==0){next}
	
	if(keep_dangling_edges)
	{
		if(n==0 && ($2 in counter) && $9>=mapq)
		{
			n=1; a[1]=$2
		}
		if(m==0 && ($6 in counter) &&$12>=mapq)
		{
			m=1; b[1]=$6
		}
	}
	
	if(report_combinatorial_edges)
	{
		delete reported 	## keep track of prev reported to track SNPs read from both sides..

		for(i=1;i<n;i++)
		{
			for(j=i+1;j<=n;j++)
			{
				$2=a[i]; $6=a[j];
				if(!keep_native_coordinates){$3=1; $7=1;};
				print; if(amplify_neighboring_snps){print};
				reported[a[i]" "a[j]]=1
			}
		}
		
		for(i=1;i<m;i++)
		{
			for(j=i+1;j<=m;j++)
			{
				if (reported[b[i]" "b[j]]){continue}
				$2=b[i]; $6=b[j];
				if(!keep_native_coordinates){$3=1; $7=1;};
				print; if(amplify_neighboring_snps){print};
				reported[b[i]" "b[j]]=1
			}
		}
				
		for(i=1;i<=n;i++)
		{
			for(j=1;j<=m;j++)
			{
				if (a[i]==b[j] || reported[a[i]" "b[j]] || reported[b[j]" "a[i]]){continue}
				$2=a[i]; $6=b[j]; 
				if(!keep_native_coordinates){$3=1; $7=1;};
				print;
			}
		}
	}
	else
	{
		if(n==0||m==0){next}

		# just take one random overlapping SNP per one read half
		n=int(1+rand()*n)
		m=int(1+rand()*m)
		$2=a[n]; $6=b[m];
		if(!keep_native_coordinates){$3=1; $7=1;};
		print
	}
}
