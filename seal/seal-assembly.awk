## Seal script to dump sealed assembly.
## awk -v size=reliable_decision_size -v bundle=bundle_size -f seal-assembly.awk <path-to-assembly-file>
## Output: assembly-formatted stdout
## Written by: OD
## Version: 180922 
## NOTE: relies on standard notation :::fragment_ and :::debris
## This version adds adds back :::debris elements that are possibly interspersed with small elements likely filtered due to scaffolder input size limit.
## NOTE: the script also sorts small non-fragment singletons and adds them to the very end to facilitate bundling, max size to bundle is passed with -v bundle option

## utility function
function delete_last(input_str,			out_str, a, n, i){
	n=split(input_str, a)
	for(i=1;i<=n-1;i++)
	{
		out_str=out_str" "a[i]
	}
	return substr(out_str, 2)
}

# read in assembly file properties section and list possible false-positives
# trying to catch Fn-Dn+1-Fn+2 and Fn-(Dn+1-Sn+2-...Sn+k-1-Dn+k)-Fn+k+1 where D stands for debris fragments and S stands for small fragments
$1~/^>/{
	split($1, a, ":::fragment_||:::debris")
	if ($1~/:::fragment_/ && ($3<size || $1~/:::debris/))
	{
		# will attempt to join things that were annotated as debris or filtered out because of size constraint
		candidate[$2]=a[1]
		candidate[-$2]=a[1]	
	}
	frag[$2]=a[1]
	frag[-$2]=a[1]
	name[$2]=$1
	len[$2]=$3
	next
}
# read in assembly file asm section
{
	# store original lines
	original_line_count++
	line[original_line_count]=$0
	
#	this is not necessary but we can keep track of unattempted to keep them at the end
	if (NF==1 && $1>0 && name[$1]!~/:::fragment_/ && len[$1]<bundle)
	{
		unattempted[$1]=$1
		next
	}
		
	# split elements consisting entirely from candidates and remove singleton orientation for uniform downstream treatment

	for(i=1; i<=NF; i++){
		if (!($i in candidate))
		{
			next
		}
	}
	gsub("-","")
	original_line_count--
	for(i=1;i<=NF;i++){
		original_line_count++
		line[original_line_count]=$i
		singleton[$i]=$i
	}
}
END{
# 1) merge candidates together

	n=asort(singleton,sorted_singleton,"@val_num_asc")
	singleton_counter=1
	str[singleton_counter]=" "sorted_singleton[1]
	for(i=2;i<=n;i++)
	{
		if (sorted_singleton[i]==sorted_singleton[i-1]+1)
		{
			if (frag[sorted_singleton[i]]==frag[sorted_singleton[i-1]])
			{
				str[singleton_counter]=str[singleton_counter]" "sorted_singleton[i]
				continue
			}
		}		
		
		str[singleton_counter]=substr(str[singleton_counter],2)
		singleton_counter++
		str[singleton_counter]=" "sorted_singleton[i]
	}
	str[singleton_counter]=substr(str[singleton_counter],2)

#  for (i in str) {print str[i] > "/dev/stderr"}
#  exit

#	2) insert merged candidates into assembly whenever possible

	new_line_count=0

	for (i=1; i<=original_line_count; i++)
	{

		if (line[i] in singleton)	# handle singletons separately
		{
			continue
		}
		
		if (line[i] in unattempted)	# handle unattempted separately, not necessary but we can keep order as in original just in case
		{
#			unattemptedstr=unattemptedstr" "line[i]
			continue
		}
		
		n=split(line[i], a)
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
						for(tmp=1;tmp<=m;tmp++){len[a[k]]+=len[b[tmp]]; delete len[b[tmp]]}; len[a[k]]+=len[a[k-1]]; delete len[a[k-1]];
						name[a[k]]=name[a[k-1]]
						output=delete_last(output)																		
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
			if(output=="")	#spaces are harmless but annoying
			{
				output=a[k]
			}
			else
			{
				output=output" "a[k]
			}
		}
		new_line_count++
		newasm[new_line_count]=output
	}
	
#	3) print remaining singletons, keep pretty order
	
	for (i=1; i<=singleton_counter; i++)
	{
		if (!str[i]) continue; # skip deleted 
		n=split(str[i], a)
		if (n>1)
		{
			for(tmp=2;tmp<=n;tmp++)
			{
				len[a[1]]+=len[a[tmp]];
				delete len[a[tmp]]
			}
			if (name[a[1]]!~/:::debris/) # add debris annotation for complex singletons
			{
				name[a[1]]=name[a[1]]":::debris"
			}
		}
		new_line_count++
		newasm[new_line_count]=a[1]
		
#		print name[newasm[new_line_count]]" "newasm[new_line_count]" "len[newasm[new_line_count]] > "/dev/stderr"
	}
#	exit
	
#	4) print remaining unattempted in sorted form - this is not necessary but useful for bundling

	#n=split(unattemptedstr,a," ")
	n=asort(unattempted)
	for (i=1;i<=n;i++)
	{	
		new_line_count++
		newasm[new_line_count]=unattempted[i]
	}
	
#	5) update properties	

	counter=0
	for (i=1; i<=length(name); i++)
	{
		if (!(i in len)){continue} # skip deleted
		
		counter++
		newid[i]=counter;
		newid[-i]=-counter;
		newname[i]=name[i]
		n=split(name[i], b, ":::fragment_||:::debris")

		# shift fragment number
		if (n>1 && name[i]~/:::fragment_/)
		{
			if (b[2]==1)
			{
				if (remember){sub(/:::fragment_1/,"", remember); print remember; remember="";}
				fragcounter=1
				remember=newname[i]" "newid[i]" "len[i]
				continue
			} else
			{
				fragcounter++
				newname[i]=b[1]":::fragment_"fragcounter
				if(name[i]~/:::debris/){newname[i]=newname[i]":::debris"}
			}
		}
				
 		if (remember){
 			if (fragcounter==2)
 			{
 				print remember; remember=""
 			} else
 			{
 				sub(/:::fragment_1/,"", remember);
 				print remember; remember="";
 			}
 		}
		print newname[i], newid[i], len[i]
	}
	# in case last line was merged
	if (remember){
 			if (fragcounter==2)
 			{
 				print remember; remember=""
 			} else
 			{
 				sub(/:::fragment_1/,"", remember);
 				print remember; remember="";
 			}
 	}
#	6) update assembly

	for (i=1; i<=new_line_count; i++)
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