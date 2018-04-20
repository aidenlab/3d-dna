## Extract highest stanza from LASTZ output.
## Written by: OD
## Version 170217.
$0=="s {"{stanza="seq"; next}
$0=="h {"{stanza="head"; next}
$0=="a {"{stanza="align"; next}
$0=="}"{stanza=""; next}
stanza!=""{

	if (stanza=="seq")
	{
		if ($(NF-1)==1)
			rev=1
		else
			rev=0
		if (len1==0)
		{
			len1=$3
			next
		}
		if (len2==0)
			len2=$3
	}
	
	if (stanza=="head")
	{
		gsub(">||\"","")
		if (!name1)
			name1=$1
		else if (!name2)
			name2=$1
	}
	
	if (stanza=="align")
	{
		if ($1=="s" && $2>=max && rev==0)
		{
			tentmax=$2
			tentative=1
		}
		if ($1=="b" && tentative==1 && $2-$3>=last_tiled_pos)
		{
			target_start=$2
			query_start=$3
			max=tentmax
		}
		if ($1=="e" && tentative==1 && max==tentmax)
		{
			target_end=$2
			query_end=$3
			tentative=0
		}
	}
}
END{
		if (max)
		{
			print target_start, target_end, len1, query_start, query_end, len2, max_score
		}
		else
		{
			print ""
			print ""
		}
		
#		if (sum[0]>sum[1])
#			rev=0
#		else
#		{
#			rev=1
#			tmp1=len2-query_start[rev]+1
#			tmp2=len2-query_end[rev]+1
#			query_start[rev]=tmp2
#			query_end[rev]=tmp1
#		}
#		print sum[rev], len1, len2, name1, target_start[rev], target_end[rev], name2, query_start[rev], query_end[rev], rev
}
