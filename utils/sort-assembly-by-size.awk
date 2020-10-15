## OD: 200616
## helper script to sort the first n sequences in the assembly by size.
## awk [-v n=seq_to_sort] -f sort-assembly-by-size.awk <path_to_asembly_file>
# read in cprops part
$0~/^>/{
	print
	len[$2]=$3
	next
}
# read in asm part
{
	counter++
	str[counter]=$0;
	c=0;
	if(!n||counter<=n){	
		gsub("-","")
		for(i=1; i<=NF; i++)
		{
			c+=len[$i]
		}
	}
	line[c]=line[c]" "counter
}
END{
	m=asorti(line, sortedline, "@ind_num_asc")
	for(i=m; i>=1; i--){
		k=split(substr(line[sortedline[i]],2), a, " ")
		for(s=1; s<=k; s++){
			print str[a[s]]
		}
	}
}
