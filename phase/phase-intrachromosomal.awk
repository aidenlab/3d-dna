## Main phasing script. Expects all SNPs in blocks file to live on one chromosome.
## USAGE: awk [ -v background=1 -v stringency=3 -v verbose=0 -v outfile=out.psf ] -f phase-intrachromosomal <path to in.psf> <path to edges.mnd.txt>
## INPUT: prephase set file describing available phased blocks (in trivial case 1 block = 1 SNP); dedupped file describing Hi-C contacts between different SNPs (edges.mnd.txt)
## OUPUT: output file in .psf format.
## PARAMETERS: background (default: 1), stringency (default: 3), verbose (default 0, i.e. false). 
## Written by O.D. (olga.dudchenko@bcm.edu)
## Version: 5/22/19
## NOTE: unprompted option ignore_sort [=0 by default i.e. false].
## NOTE: passing n as argument for verbose will report progress every nth step.

function invert(string,		str,i,n)
{
	if (string=="")
		return ""
	n=split(string, entries, " ")
	for(i=1;i<=n;i++)
	{
		str=str" "(-entries[i])
	}
	return substr(str,2)
}
BEGIN{
	# set defaults
	if(!background){background = 1}
	if(!stringency){stringency = 1/3}else{stringency = 1/stringency}
	if(!outfile){outfile="out.psf"}

	megacyclecounter=0
}
## read in .psf file, properties section: this describes individual SNPs. Assumed format: ">CHR POS REF ALT ID"
FILENAME==ARGV[1]&&($1~/^>/){
	
	properties[NR]=$0
	participating[$NF]=1
	if (length(chr)==0){chr=$1}else{if($1!=chr){print ":( SNPs span multiple chromosomes. Exiting!" > "/dev/stderr"; force_exit = 1; exit}}
	
	$1=substr($1,2)
		
	snp[$1":"$2":"$3]=$NF
	snp[$1":"$2":"$4]=$NF
	phase[$1":"$2":"$3]=1
	phase[$1":"$2":"$4]=-1
	
	next
}
## read in .psf file, set section: this describes prephased blocks. In trivial case 1 SNP = 1 block.
FILENAME==ARGV[1]{
	counter++
	split($0, a, " ")
	for (i in a)
	{
		if ((!(a[i] in participating)) && (!(-a[i] in participating))){print ":( Blocks list SNPs not described in the properties section. Exiting!" > "/dev/stderr"; force_exit = 1; exit}
		if (a[i]>0)
		{
			blockId[a[i]]=counter
			inBlockPhase[a[i]]=1
		}
		else
		{
			blockId[-a[i]]=counter
			inBlockPhase[-a[i]]=-1		
		}
	}
	block[counter]=$0	# trivial block[n]=n
	next
}
## read in pairwise SNP edge file, count pscores and nscores, for individual snp pairs
(!($2 in snp)) || (!($6 in snp)){next}
{
	if (phase[$2]*phase[$6] > 0)
		pscore[snp[$2]" "snp[$6]]++
	else
		nscore[snp[$2]" "snp[$6]]++
}
## actual phasing
END{
	
## handle exit	
	if (force_exit){exit 1}

	print ":) Phasing chromosome "substr(chr,2)"." > "/dev/stdout"	

## INITIALIZATION

## from pairwise SNP scores to pairwise block scores, collect partners and calculate first set of confidences: catch up with prephased blocks or trivial pscore=pbscore; nscore=nbscore

	print "... :) Initializing blocks." > "/dev/stdout"

	### from snps to blocks
	for (i in pscore)
	{
		split(i, a, " ")
		if (blockId[a[1]]==blockId[a[2]])
			continue
		else if (blockId[a[1]] < blockId[a[2]])
			var=blockId[a[1]]" "blockId[a[2]]
		else
			var=blockId[a[2]]" "blockId[a[1]]
			
		if (! (var in pbscore))
		{
			partners[blockId[a[1]]]=partners[blockId[a[1]]]" "blockId[a[2]]
			partners[blockId[a[2]]]=partners[blockId[a[2]]]" "blockId[a[1]]
		}
		
		if (inBlockPhase[a[1]]==inBlockPhase[a[2]])
		{
			pbscore[var]+=pscore[a[1]" "a[2]]
			nbscore[var]+=0
		}
		else
		{
			nbscore[var]+=pscore[a[1]" "a[2]]
			pbscore[var]+=0
		}

	}
	
	for (i in nscore)
	{
		split(i, a, " ")
		if (blockId[a[1]]==blockId[a[2]])
			continue
		else if (blockId[a[1]] < blockId[a[2]])
			var=blockId[a[1]]" "blockId[a[2]]
		else
			var=blockId[a[2]]" "blockId[a[1]]

		if(!(var in pbscore))
		{
			partners[blockId[a[1]]]=partners[blockId[a[1]]]" "blockId[a[2]]
			partners[blockId[a[2]]]=partners[blockId[a[2]]]" "blockId[a[1]]
		}
			
		if (inBlockPhase[a[1]]==inBlockPhase[a[2]])
		{
			nbscore[var]+=nscore[a[1]" "a[2]]
			pbscore[var]+=0
		}
		else
		{
			pbscore[var]+=nscore[a[1]" "a[2]]
			nbscore[var]+=0
		}
			
	}
	
	### add trailing space in partners, to play a trick with substitution; format _1_2_3_
	for (i in partners)
	{
		partners[i]=partners[i]" "
	}

	### from now on, everything should be in terms of blocks, not individual snps
	
	### confidences for all block pairs
	for (i in pbscore)	# pbscore and nbscore defined for the same set, does not matter
	{
		diff=pbscore[i]-nbscore[i]
		codirected[i]=1
		if(diff<0)
		{
			diff=-diff
			codirected[i]=0
		}
		sum[i]=pbscore[i]+nbscore[i]
		tmp=(diff/stringency/(sum[i]+background))
		if (tmp>1)
			conf[i]=tmp	
	}

## MAIN CYCLE

	print "... :) Starting main program cycle. Block count coming in: "length(block)"." > "/dev/stdout"
	
	## merge blocks and update appropriate variables, starting with pairs of highest confidence
	while(length(conf)>0)
	{	
	
		if(verbose && !(megacyclecounter%verbose))
		{
			### print intermediate results	
			for(i=1;i<=length(properties); i++)
				print properties[i] > "verbose."ARGV[1]"."megacyclecounter".psf"
			### sort by block number and within blocks, by default...
			for(i=1;i<=length(properties); i++)
			{
				if(block[i])
				{
					blockcounter++
					if (ignore_sort) # i don't know why i am keeping this
						print block[i] > "verbose."ARGV[1]"."megacyclecounter".psf"
					else
					{	
						print block[i] > "verbose."ARGV[1]".tempfile"
						command="cat verbose."ARGV[1]".tempfile | tr \" \" \"\\n\" | awk '$1>0{print $1, $1;next}{print -$1, $1}' | sort -k 1,1n | awk '{print $2}' | paste -s -d \" \" >> verbose."ARGV[1]"."megacyclecounter".psf"
						system(command)
						close(command)
						close("verbose."ARGV[1]".tempfile")
					}
				}
			}

# 			for(i in block)
# 			{
# 				print block[i] > "h."megacyclecounter".psf"
# 			}
		print "	...iterative cycle #"megacyclecounter". Block count in this iteration: "blockcounter"." > "/dev/stdout"
		}
		
		megacyclecounter++
		# if (verbose && !(megacyclecounter%verbose))
		# 	print "	...iterative cycle #"megacyclecounter". Block count in this iteration: "blockcounter"." > "/dev/stdout"

	
# 0) clean temporary arrays
		split("", modified)
		split("", sorted_pair)
 
	
# 1) sort block pairs by confidences and, if same, by scores. Tried within awk but much faster to do in shell.
		
		for (i in conf)
		{
			print i, conf[i], sum[i] > ARGV[1]".tempfile"
		}
		command = "sort -k 3,3nr -k 4,4nr "ARGV[1]".tempfile > "ARGV[1]".tempsortedfile"
		system(command)
		n=0
		name=ARGV[1]".tempsortedfile"
		while ((getline var < name) > 0)
		{
			n++
			split(var, a, " ")
            sorted_pair[n]=a[1]" "a[2]
        }
        close(command)
        close(ARGV[1]".tempfile")
        close(ARGV[1]".tempsortedfile")
			
# 2) merge pairs of blocks, prioritizing according to sorting from prev step				
		
		for(i=1; i<=n; i++)
		{				
			if(!(sorted_pair[i] in conf)){continue}	# e.g. second in pair was modified and associated confidences already merged
													
			split(sorted_pair[i], a, " ")	# get block numbers - merge candidates
					
			### if block has been already updated in this cycle, skip. Note this is a conservative variation of this script: keep track of things that have a chance of being updated with higher confidence...
			if ((a[1] in modified) || (a[2] in modified)){modified[a[1]]=1; modified[a[2]]=1; continue} # a[2] should be handled by length check #228

			### actually go for a merge

			modified[a[1]]=1; modified[a[2]]=1; # keeping track of a[2] is excessive, clean up later?
			
			### delete confidence of the pair to be merged
			delete conf[sorted_pair[i]]	# remove confidence associated with current block pair
			
			### remove the second block from list of partners
			gsub(" "a[2]" ", " ", partners[a[1]])

			if(codirected[sorted_pair[i]]==1)	### merge codirected
			{
				# update block			
				block[a[1]]=block[a[1]]" "block[a[2]]
				delete block[a[2]]

				# iterate through all the partners of the second block
				split(substr(partners[a[2]],2, length(partners[a[2]])-2), secondblockpartners, " ")
				for (k in secondblockpartners)
				{
					if (secondblockpartners[k]==a[1])
						continue						
					
					if (a[1]<secondblockpartners[k])
						tmpname=a[1]" "secondblockpartners[k]
					else
						tmpname=secondblockpartners[k]" "a[1]
					
					# update partner lists				
					if (!(tmpname in pbscore))
					{
						partners[a[1]]=partners[a[1]]""secondblockpartners[k]" "
						gsub(" "a[2]" ", " "a[1]" ", partners[secondblockpartners[k]])						
					}
					else
					{
						gsub(" "a[2]" " ," ", partners[secondblockpartners[k]])
					}
						
					# update scores
					if (a[2]<secondblockpartners[k])
					{												
						pbscore[tmpname]+=pbscore[a[2]" "secondblockpartners[k]]
						nbscore[tmpname]+=nbscore[a[2]" "secondblockpartners[k]]
						delete pbscore[a[2]" "secondblockpartners[k]]
						delete nbscore[a[2]" "secondblockpartners[k]]
						delete conf[a[2]" "secondblockpartners[k]]	# if existed
					}
					else
					{							
						pbscore[tmpname]+=pbscore[secondblockpartners[k]" "a[2]]
						nbscore[tmpname]+=nbscore[secondblockpartners[k]" "a[2]]							
						delete pbscore[secondblockpartners[k]" "a[2]]
						delete nbscore[secondblockpartners[k]" "a[2]]								
						delete conf[secondblockpartners[k]" "a[2]]	# if existed
					}

					# update confidences of remaining blocks
					delete conf[tmpname]	# if existed
				
					diff=pbscore[tmpname]-nbscore[tmpname]
					codirected[tmpname]=1
					if(diff<0)
					{
						diff=-diff
						codirected[tmpname]=0
					}
					sum[tmpname]=pbscore[tmpname]+nbscore[tmpname]
					tmp=diff/stringency/(sum[tmpname]+background) # something happening here when switching to >1 stringency, need to explore
					if (tmp>1)
						conf[tmpname]=tmp
				}
				delete partners[a[2]]
			}
			else	### merge antidirected
			{
				# update block			
				block[a[1]]=block[a[1]]" "invert(block[a[2]])
				delete block[a[2]]

				# iterate through all the partners of the second block
				split(substr(partners[a[2]],2, length(partners[a[2]])-2), secondblockpartners, " ")
				for (k in secondblockpartners)
				{
					if (secondblockpartners[k]==a[1])
						continue						
					
					if (a[1]<secondblockpartners[k])
						tmpname=a[1]" "secondblockpartners[k]
					else
						tmpname=secondblockpartners[k]" "a[1]
					
					# update partner lists				
					if (!(tmpname in pbscore))
					{
						partners[a[1]]=partners[a[1]]""secondblockpartners[k]" "
						gsub(" "a[2]" ", " "a[1]" ", partners[secondblockpartners[k]])						
					}
					else
					{
						gsub(" "a[2]" " ," ", partners[secondblockpartners[k]])
					}
										
					# update scores
					if (a[2]<secondblockpartners[k])
					{												
						pbscore[tmpname]+=nbscore[a[2]" "secondblockpartners[k]]
						nbscore[tmpname]+=pbscore[a[2]" "secondblockpartners[k]]
						delete pbscore[a[2]" "secondblockpartners[k]]
						delete nbscore[a[2]" "secondblockpartners[k]]
						delete conf[a[2]" "secondblockpartners[k]]	# if existed
					}
					else
					{							
						pbscore[tmpname]+=nbscore[secondblockpartners[k]" "a[2]]
						nbscore[tmpname]+=pbscore[secondblockpartners[k]" "a[2]]							
						delete pbscore[secondblockpartners[k]" "a[2]]
						delete nbscore[secondblockpartners[k]" "a[2]]								
						delete conf[secondblockpartners[k]" "a[2]]	# if existed
					}

					# update confidences of remaining blocks
					delete conf[tmpname]	# if existed
				
					diff=pbscore[tmpname]-nbscore[tmpname]
					codirected[tmpname]=1
					if(diff<0)
					{
						diff=-diff
						codirected[tmpname]=0
					}
					sum[tmpname]=pbscore[tmpname]+nbscore[tmpname]
					tmp=diff/stringency/(sum[tmpname]+background)
					if (tmp>1)
						conf[tmpname]=tmp
				}
				delete partners[a[2]]
			}
		}
	}

## PRINT RESULTS	

	for(i=1;i<=length(properties); i++)
		print properties[i] > outfile
	
	# sort by block number and within blocks.
	for(i=1;i<=length(properties); i++)
	{
		if(block[i])
		{
			if (ignore_sort) # i don't know why i am keeping this
				print block[i] > outfile
			else
			{	
				blockcounter++		
				print block[i] > ARGV[1]".tempfile"
				command="cat "ARGV[1]".tempfile | tr \" \" \"\\n\" | awk '$1>0{print $1, $1;next}{print -$1, $1}' | sort -k 1,1n | awk '{print $2}' | paste -s -d \" \" >> "outfile				
				
				
				system(command)
				close(command)
				close(ARGV[1]".tempfile")
			}
		}
	}
	system("rm -f "ARGV[1]".tempfile "ARGV[1]".tempsortedfile")
	if (verbose)
		system("rm -f verbose."ARGV[1]".tempfile verbose."ARGV[1]".tempsortedfile verbose."ARGV[1]".*.psf")
	print "... :) Done with main program cycle! Block count coming out: "blockcounter"." > "/dev/stdout"
}
