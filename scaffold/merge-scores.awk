#### Description: Normalize raw counts for density graph. Collect data in case input was parallelized.
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
		}
		next
	}
}
# read in the sorted score stdin
{
	if (($1!=prev1 || $2!=prev2 || $3!=prev3) && FNR!=1)
	{
		print prev1, prev2, prev3, count/(slength[prev1]*slength[prev2]), count
		count = 0
	}
	prev1=$1; prev2=$2; prev3=$3
	count += $4
}
END{
	if (prev1)
	{
		print prev1, prev2, prev3, count/(slength[prev1]*slength[prev2]), count
	}
}