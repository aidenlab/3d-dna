## OD: 171118
# read in cprops
FILENAME==ARGV[1]{
	if ($1!~/:::fragment_/ && $3<input_size){$3=0}
	len[$2]=$3
	next
}
# read in asm
{
	str[FNR]=$0;
	c=0;
	gsub("-","")
	for(i=1; i<=NF; i++){
		c+=len[$i]
	}
	line[c]=line[c]" "FNR
}
END{
	n=asorti(line, sortedline, "@ind_num_asc")
	for(i=n; i>=1; i--){
		k=split(substr(line[sortedline[i]],2), a, " ")
		for(s=1; s<=k; s++){
			print str[a[s]]
		}
	}
}