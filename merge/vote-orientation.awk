## Part of the merge segment of the diploid pipeline. Decides the relative orientation of the connected component with respect to the rest of the assembly by majority vote.
## Written by: OD
## read in cprops
FILENAME==ARGV[1]{
	cname[$1]=$2
	oname[$2]=$1
	len[$1]=$3
	next
}
## read in asm
FILENAME==ARGV[2]{
	for(i=1;i<=NF;i++)
	{
		if($i>0)
			orientation[$i]=1
		else
			orientation[-$i]=-1
	}
	next
}
## read cluster
{
	len_total=len[substr($1,2)]
	len_conflicting=0
	delete conflicting
	for(i=2;i<=NF;i++)
	{
		len_total+=len[substr($i,2)]
		if ((substr($i,1,1)=="+" && orientation[cname[substr($i,2)]]!=orientation[cname[substr($1,2)]]) || (substr($i,1,1)=="-" && orientation[cname[substr($i,2)]]==orientation[cname[substr($1,2)]]))
		{
			len_conflicting+=len[substr($i,2)]
			conflicting[$i]=1
		}
	}	
	
	if (len_conflicting <= len_total/2)
	{
		str=orientation[cname[substr($1,2)]]*cname[substr($1,2)]
		for(i=2;i<=NF;i++)
		{
			if ($i in conflicting)
				str=str" "(-1)*orientation[cname[substr($i,2)]]*cname[substr($i,2)]
			else
				str=str" "orientation[cname[substr($i,2)]]*cname[substr($i,2)]
		}
	}
	else
	{
		str=(-1)*orientation[cname[substr($1,2)]]*cname[substr($1,2)]
		for(i=2;i<=NF;i++)
		{
			if ($i in conflicting)
			{
				str=str" "orientation[cname[substr($i,2)]]*cname[substr($i,2)]
			}
			else
			{
				str=str" "(-1)*orientation[cname[substr($i,2)]]*cname[substr($i,2)]
			}
		}
	}
	print str
}