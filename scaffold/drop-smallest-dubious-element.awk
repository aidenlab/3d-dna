#### Description: helper script to remove the smallest scaffolds from the scaffold pool.
#### Written by: OD
BEGIN{
	smallest=0
}
function abs(value)
{
  return (value<0?-value:value);
}
# read in the tentative dropouts list
{
	if (FILENAME==ARGV[1])
	{
		dropout[$1]=1
		next
	}
}
# read in the cprops
{
	if (FILENAME==ARGV[2]) 
	{
		cname[$1]=$2
		clength[$2]=$3
		next
	}

}
# figure out the smallest dropout candidate
{
	if (FNR in dropout) 
	{
		n=split($0,a)
		slength[FNR]=0
		for (k=1; k<=n; k++)
		{
			slength[FNR]+=clength[abs(a[k])]
		}
		if ((smallest==0)||(slength[FNR]<=smallest))
		{
			smallest=slength[FNR]
			candidate=FNR
		}
		
	}
}
END{
	print candidate
   }

