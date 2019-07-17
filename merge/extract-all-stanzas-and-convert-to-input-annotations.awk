## Extract highest stanza from LASTZ output. Only one orientation is allowed.
## Options: skip_header={0,1}
## Written by: OD

BEGIN{
	OFS="\t"
	if(!skip_header){print "chr1	x1	x2	chr2	y1	y2	color	id	X1	X2	Y1	Y2"}
}

$0=="s {"{stanza="seq"; len1=0; len2=0; next}
$0=="h {"{stanza="head"; target_name=""; query_name=""; next}
$0=="a {"{stanza="align"; next}
$0=="}"{stanza=""; next}
stanza!=""{

	if (stanza=="seq")
	{
		if ($(NF-1)==1)
			rev=1
		else
			rev=0
		if (!len1)
		{
			len1=$3
			next
		}
		if (!len2)
		{
			len2=$3
			next	
		}
	}
	
	if (stanza=="head")
	{
		gsub(">||\"","")
		if (!target_name)
		{
			target_name=$1
			next
		}
		if (!query_name)
		{
			query_name=$1
			next
		}
	}
	
	if (stanza=="align")
	{
		if ($1=="s")
		{
			#counter++
			score=$NF
			#orientation=rev
			next
		}
		if ($1=="b")
		{
			target_begin=$2
			query_begin=$3
			next
		}
		if ($1=="e")
		{
			target_end=$2
			query_end=$3
			
			# dump results:
			
			if(rev==1){
				tmp1=len2-query_begin+1
				tmp2=len2-query_end+1
				query_begin=tmp2
				query_end=tmp1
				score="-"score
			}else{
				score="+"score
			}
		}
		
		print target_name, target_begin, target_end, query_name, query_begin, query_end, "0,0,0", score, target_begin, target_end, query_begin, query_end
		# in case this is a combined file
		stanza=""
	}
}
