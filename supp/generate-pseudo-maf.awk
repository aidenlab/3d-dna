## Helper script to generate maf-looking files for before/after assembly dotplots.
## awk -c <chra|chrb> -r <100000> -f generate-pseudo-maf.awk <before.assembly>
## OD: version 210309
BEGIN{
	if(length(c)){split(c,a,"|"); for(i in a){list[a[i]]=1}}
	if(!length(r)){res=100000}else{res=r}
}
$0~/^>/{chr=substr($1,2); if(c&&(!(chr in list))){next}; list[chr]=1; len[chr]=$3}
END{
	for(i in list)
	{
		if (len[i]<res){continue}
		pos=1;
		while(pos <= len[i]-res)
		{
			print "a score="res
			print "s "i" "pos" "res" + "len[i]" N"
			print "s "i" "pos" "res" + "len[i]" N"
			print ""
			pos+=res
		}
	}

}
