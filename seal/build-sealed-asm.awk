## Seal script to dump new asm. Creates lits of annotations to be deleted (filename1) and to be merged (filename2) in the process
## awk -v <filename1> -v <filename2> -f build-sealed-asm.awk <path-to-cprops> <path-to-asm>
## Output: asm-formatted stdout and some additional temporary files
## Written by: OD
## NOTE: relies on standard notation
function delete_last(input_str,			out_str, a, n, i){
	n=split(input_str, a)
	for(i=1;i<=n-1;i++)
	{
		out_str=out_str" "a[i]
	}
	return substr(out_str, 2)
}
BEGIN{
	if (! filename1)
		filename1="h.to_delete.txt"
	if (! filename2)
		filename2="h.to_merge.txt"
}
# read in cprops
FILENAME==ARGV[1]{
	split($1, a, ":::fragment_||:::debris")
	if ($1~/:::fragment_/ && ($3<size || $1~/:::debris/))
	{
		candidate[$2]=a[1]	# will attempt to join things that were annotated as debris or filtered out because of size constraint. In principle could have join all singletons, dropouts are rare...
	}
	frag[$2]=a[1]
	frag[-$2]=a[1]
	next
}
# read in asm
{
	line[FNR]=$0
	if (NF==1)
		{
			gsub("-","")
			singleton[$1]=$1
		}
}
END{

	# list connected candidates
	n=asort(singleton)
	counter=1
	str[counter]=" "singleton[1]
	for(i=2;i<=n;i++)
	{
		if (singleton[i]==singleton[i-1]+1)
		{
			if ((singleton[i] in candidate) && (singleton[i-1] in candidate) && (candidate[singleton[i]]==candidate[singleton[i-1]]))
			{
				str[counter]=str[counter]" "singleton[i]
				continue
			}
		}		
		
		str[counter]=substr(str[counter],2)
		counter++
		str[counter]=" "singleton[i]
	}
	str[counter]=substr(str[counter],2)

	# build new asm in terms of old cprops

	line_count=0

	for (i=1; i<=FNR; i++)
	{
		n=split(line[i], a)
		if (n==1)	# handle singletons separately
		{
			continue		
		}
		
		output=a[1]
		
		for(k=2; k<=n; k++)
		{

			if (frag[a[k]]==frag[a[k-1]])
			{
				for (s in str)
				{
					m=split(str[s],b)
					if (a[k-1]==b[1]-1 && a[k]==b[m]+1)
					{
						output=delete_last(output)						
						a[k]=a[k-1]
						print str[s] > filename1
						delete str[s]
						break
					}
					else if (a[k-1]==-(b[m]+1) && a[k]==-(b[1]-1))
					{
						output=delete_last(output)
						print str[s] > filename1
						delete str[s]
						break
					}
				}
			}
			output=output" "a[k]		
		}
		line_count++
		newasm[line_count]=output
	}
	# print remaining singletons
	for (i in str)
	{
		n=split(str[i], a)
		if (n>1)
		{
			print str[i] > filename2
		}
		
		line_count++
		newasm[line_count]=a[1]
	}	
	
	## shift numbers
	for (i=1; i<=line_count; i++)
	{
		all_scaffolds=all_scaffolds" "newasm[i]
	}
	gsub("-","",all_scaffolds)
	n=split(all_scaffolds, a)
	asort(a)
	for(i=1; i<=length(a); i++)
	{
		newid[a[i]]=i
		newid[-a[i]]=-i
	}
	
	## print asm
	for (i=1; i<=line_count; i++)
	{
		n=split(newasm[i],a)
		out_str=""
		for(k=1; k<=n; k++)
		{
			out_str=out_str" "newid[a[k]]
		}
		print substr(out_str,2)
	}
}