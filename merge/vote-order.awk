## Part of the merge segment of the diploid pipeline. Decides the relative position of the connected components with respect to the rest of the assembly by majority vote.
## Written by: OD
# read in cprops
FILENAME==ARGV[1]{
	len[$2]=$3
	len[-$2]=$3
	next
}
# read in clusters
FILENAME==ARGV[2]{
	cluster[$0]=$0
	next
}
# read in asm line by line
{
	delete shift
	shift[$1]=0
	for(i=2;i<=NF;i++)
	{
		shift[$i]=shift[$(i-1)]+len[$(i-1)]
	}
	for (i in cluster)
	{
		n=split(i, a)
		if (shift[a[1]]=="" && shift[-a[1]]=="")
		{
			continue
		}
		for(k=1; k<=n; k++)
		{
			if (shift[a[k]])
			{
				shift[i]+=len[a[k]]*shift[a[k]]
			}
			else
			{
				shift[i]+=len[a[k]]*shift[-a[k]]
			}
			delete shift[a[k]]
			delete shift[-a[k]]
			len[i]+=len[a[k]]
		}
		shift[i]=int(shift[i]/len[i])
	}
	
	n=asort(shift, dest)	
	outstr=""
	for(i=1; i<=n; i++)
	{
		for(k in shift)
		{
			if(shift[k]=="")
				continue
			if(dest[i]==shift[k])
			{
				if (k in cluster)
					outstr=outstr" {"k"}"
				else
					outstr=outstr" "k
				delete shift[k]			
			}
		}
	}
	print substr(outstr,2)
}