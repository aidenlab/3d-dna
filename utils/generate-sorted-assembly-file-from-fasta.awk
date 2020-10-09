#### Utility script to generate .assembly file from fasta. The order of scaffolds in the .assembly file is from largest to smallest.
#### USAGE: awk -f generate-sorted-assembly-file-from-fasta.awk <path-to-fasta>
#### Output: assembly-formatted stdout.

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
			print ">"a[s], counter, sortlen[i]
		}
	}
	for(i=1; i<=counter; i++){
		print i
	}
}
