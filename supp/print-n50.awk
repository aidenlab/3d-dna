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
	printf("Total contig length:\t%'d\n", totallength)
	printf("No of contigs:\t%'d\n", s)
	printf("contig N50:\t%'d\n", N50)
	printf("Largest contig:\t%'d\n", clength[s])
        
        
    print "---------"
        
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
	printf("Total scaffold length:\t%'d\n", totallength)
	printf("No of scaffolds:\t%'d\n", s)
	printf("scaffold N50:\t%'d\n", N50)
	printf("Largest scaffold:\t%'d\n", slength[s])

}
