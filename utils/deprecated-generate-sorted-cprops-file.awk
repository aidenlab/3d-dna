#### Description: Generate cprops file describing one-to-one correspondence between original draft sequences names and internal naming convention. Includes total lengths.
#### USAGE: awk -f generate-cprops-file.awk <path-to-fasta>
#### Output: cprops-formatted stdout
$0~/>/{
	if (name && c){len[c]=len[c]" "name}
	name=substr($1,2);
	c=0
	next
}
{
	c += length
}
END{
	len[c]=len[c]" "name
	n = asorti(len, sortlen, "@ind_num_asc")
	for(i=n; i>=1; i--)
	{
		k=split(substr(len[sortlen[i]],2), a, " ")
		for(s=1; s<=k; s++){
			counter++
			print a[s], counter, sortlen[i]
		}
	}
}
