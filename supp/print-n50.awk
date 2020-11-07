{
	if ($0~/>/)	# new scaffold
	{
		# contig stuff
		if (gapless_counter>0)
		{
			s+=1
			clength[s]=gapless_counter
		}
		gapless_counter=0
		
		# scaffold stuff
		if(NR>1){
			scaffold_counter++
			slength[scaffold_counter]=total_length
			total_length=0
		}
		next
	}
	
	gsub("n","N")

	total_length+=length
	
	if ($0!~/N/)
	{
		gapless_counter+=length($0)
	}
	else
	{
		n=split($0, a, "")
		for (i=1; i<=n; i++)
		{
			if (a[i]=="N")
			{
				if (gapless_counter>0)
				{
					s+=1
					clength[s]=gapless_counter
#					print s, clength[s]
				}
				gapless_counter=0
			}
			else
			{
				gapless_counter+=1
			}
		}
	}	
}
END{
    if (gapless_counter>0)
	{
		s+=1
		clength[s]=gapless_counter
	}
	
	scaffold_counter++
	slength[scaffold_counter]=total_length
	
	# print contig stuff
	asort(clength)
	for (i=1; i<=s; i++)
			totallength+=clength[i]
	i=s
	N50=clength[s]
	longerlength=clength[s]
	while (longerlength<totallength/2)
	{
			i=i-1
#               print i
			N50=clength[i]
			longerlength+=clength[i]
	}
	printf("\t\t\"totalContigLength\": \"%'d\",\n", totallength)
	printf("\t\t\"numberOfContigs\": \"%'d\",\n", s)
	printf("\t\t\"contigN50\": \"%'d\",\n", N50)
	printf("\t\t\"longestContig\": \"%'d\",\n", clength[s])
        
        
    # print scaffold stuff
    
	asort(slength)
	totallength=0
    for (i in slength)
            totallength+=slength[i]
	s=scaffold_counter
	N50=slength[s]
	longerlength=slength[s]
	i=s
	while (longerlength<totallength/2)
	{
		i--
		N50=slength[i]
		longerlength+=slength[i]
	}

        printf("\t\t\"totalScaffoldLength\": \"%'d\",\n", totallength)
        printf("\t\t\"numberOfScaffolds\": \"%'d\",\n", s)
        printf("\t\t\"scaffoldN50\": \"%'d\",\n", N50)
        printf("\t\t\"longestScaffold\": \"%'d\",\n", slength[s])
}
