#### Description: Script to build raw contact count graph.
#### Written by: OD
function abs(value)
{
	return (value<0?-value:value);
}
# read in the cprops
{
	if (FILENAME==ARGV[1]) 
	{
		cname[$1]=$2
		clength[$2]=$3
		next
	}
}
# read in current assembly
{
	if (FILENAME==ARGV[2]) 
	{
		n=split($0,a)
		slength[FNR]=0
		for (k=1; k<=n; k++)
		{
			slist[abs(a[k])]=FNR
			shift[abs(a[k])]=slength[FNR]
			slength[FNR]+=clength[abs(a[k])]
			if (a[k]==abs(a[k]))
				orientation[a[k]]=0
			else
				orientation[abs(a[k])]=1
		}
		next
	}
}
{
	#check that both reads don't fall on the same contig and have been aligned with MAPQ>=MAPQ variable - optional
	if ((cname[$2] in slist) && (cname[$6] in slist) && (slist[cname[$2]]!=slist[cname[$6]]) && ($9 >= MAPQ) && ($12 >= MAPQ)) 
	{
		$2=cname[$2]
		$6=cname[$6]			
		#recalculate positions
		if (orientation[$2] == 0)
			pos1 = $3 + shift[$2]
		else
			pos1 = -$3 + shift[$2] + clength[$2] + 1
		if (orientation[$6] == 0)
			pos2 = $7 + shift[$6]
		else
			pos2 = -$7 + shift[$6] + clength[$6] + 1
	
		## ignore everything that is not fixed length head or tail
		if ( pos1 > SIZE/2 && pos1 < slength[slist[$2]] - SIZE/2 )
			next
		if ( pos2 > SIZE/2 && pos2 < slength[slist[$6]] - SIZE/2)
			next
		
		## otherwise classify contact
		test1 = (2*pos1 - 1 - slength[slist[$2]])
		test2 = (2*pos2 - 1 - slength[slist[$6]])

		if ( test1 == 0 || test2 == 0 )
			next

		if (slist[$2] < slist[$6])
		{
			if ( test1 < 0 && test2 < 0 )
				count[slist[$2]" "slist[$6]" "1]+=1
			else if ( test1 < 0 && test2 > 0 )
				count[slist[$2]" "slist[$6]" "2]+=1
			else if ( test1 > 0 && test2 < 0 )
				count[slist[$2]" "slist[$6]" "3]+=1
			else
				count[slist[$2]" "slist[$6]" "4]+=1
		}
		else
		{
			if ( test1 < 0 && test2 < 0 )
				count[slist[$6]" "slist[$2]" "1]+=1
			else if ( test1 < 0 && test2 > 0 )
				count[slist[$6]" "slist[$2]" "3]+=1
			else if ( test1 > 0 && test2 < 0 )
				count[slist[$6]" "slist[$2]" "2]+=1
			else
				count[slist[$6]" "slist[$2]" "4]+=1
		}
	}

}
END{
for (var in count)
	{
		print var, count[var]
	}
}
