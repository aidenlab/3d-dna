## Seal script to dump sealed asm and cprops.
## awk -f build-sealed-asm-new.awk <path-to-current-cprops> <path-to-current-asm>
## Output: asm-formatted stdout and cprops-formatted stderr
## Written by: OD
## NOTE: relies on standard notation :::fragment_ and :::debris

function delete_last(input_str,			out_str, a, n, i){
	n=split(input_str, a)
	for(i=1;i<=n-1;i++)
	{
		out_str=out_str" "a[i]
	}
	return substr(out_str, 2)
}
# read in cprops
FILENAME==ARGV[1]{
	split($1, a, ":::fragment_||:::debris")
	if ($1~/:::fragment_/ && ($3<size || $1~/:::debris/))
	{
		candidate[$2]=a[1]
		candidate[-$2]=a[1]	# will attempt to join things that were annotated as debris or filtered out because of size constraint. In principle could have join all singletons, dropouts are rare...
	}
	frag[$2]=a[1]
	frag[-$2]=a[1]
	name[$2]=$1
	len[$2]=$3
	next
}
# read in asm
{
	line[FNR]=$0
	if (NF==1 && ($1 in candidate))	# leave alone big singletons
		{
			gsub("-","")
			singleton[$1]=$1
		}
}
END{
	# list connected candidates
	n=asort(singleton)
	
	# rejoin singletons
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
		if (n==1 && (a[1] in candidate))	# handle singletons separately
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
						
						for(tmp=1;tmp<=m;tmp++){len[a[k-1]]+=len[b[tmp]]; delete len[b[tmp]]}; len[a[k-1]]+=len[a[k]]; delete len[a[k]];

						output=delete_last(output)						
						a[k]=a[k-1]						
																		
						delete str[s]
						break
					}
					else if (a[k-1]==-(b[m]+1) && a[k]==-(b[1]-1))
					{
						for(tmp=1;tmp<=m;tmp++){len[-a[k]]+=len[b[tmp]]; delete len[b[tmp]]}; len[-a[k]]+=len[-a[k-1]]; delete len[-a[k-1]];

						output=delete_last(output)						
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
			for(tmp=2;tmp<=n;tmp++){len[a[1]]+=len[a[tmp]]; if(name[a[tmp]]~/:::debris/&&name[a[1]]!~/:::debris/){name[a[1]]=name[a[1]]":::debris"}; delete len[a[tmp]]}
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
	
	for(i=1;i<=length(a);i++){
		newname[a[i]]=name[a[i]]
		
		n=split(name[a[i]], b, ":::fragment_||:::debris")

 		if (n>1 && name[a[i]]~/:::fragment_/)
 		{
 			if (b[2]==1)
 			{
 				if (remember){sub(/:::fragment_1/,"", remember); print remember>"/dev/stderr"; remember="";}
 				fragcounter=1
 				remember=newname[a[i]]" "newid[a[i]]" "len[a[i]]
 				continue
 			} else
 			{
 				fragcounter++
 				newname[a[i]]=b[1]":::fragment_"fragcounter
 				if(name[a[i]]~/:::debris/){newname[a[i]]=newname[a[i]]":::debris"}
 			}
 		}
 		if (remember){
 			if (fragcounter==2)
 			{
 				print remember > "/dev/stderr"; remember=""
 			} else
 			{
 				sub(/:::fragment_1/,"", remember);
 				print remember>"/dev/stderr"; remember="";
 			}
 		}
		print newname[a[i]], newid[a[i]], len[a[i]] > "/dev/stderr"

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