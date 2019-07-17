## Parse alignment information to group contigs/scaffolds into clusters
## Written by: OD
## NOTE: if alignments are in conflict both orientations will be listed in the cluster! I.e. two entries per cluster
function add_reverse(str,	newstr,n,tmp,i){
	n=split(str,tmp)
	newstr="-"substr(tmp[1],2)
	for(i=2;i<=n;i++)
	{
		if (substr(tmp[i],1,1)=="+")
			newstr=newstr" -"substr(tmp[i],2)
		else
			newstr=newstr" +"substr(tmp[i],2)
	}
	return newstr
}
function merge(vertex,	cluster, tmp, n, i, k, c)
{
	n=split(neigh[vertex], tmp)
	cluster=tmp[1]
	for (i=2; i<=n; i++)
	{
		if (neigh[substr(tmp[i],2)] && neigh[substr(tmp[i],2)]!="")
		{
			if (substr(tmp[i],1,1)=="-")
			{
				cluster=cluster" "add_reverse(neigh[substr(tmp[i],2)])
			}
			else
			{
				cluster=cluster" "neigh[substr(tmp[i],2)]
			}
			delete neigh[substr(tmp[i],2)]
		}
		else
			cluster=cluster" "tmp[i]
	}
	# delete duplicates from cluster
	n=split(cluster, tmp)
	cluster=""
	delete c
	for(i=1;i<=n;i++)
	{
		c[tmp[i]]++
		if (c[tmp[i]]==1)
			cluster=cluster" "tmp[i]
	}
	cluster=substr(cluster,2)
	
#print "old: "neigh[vertex]
#print "new: "cluster
	if (cluster!=neigh[vertex])
	{
#print "growing..." 
		neigh[vertex]=cluster
		merge(vertex)
	}
}
## Read in the annotation file
NR>1{
	if (! neigh[$1])
	{
		neigh[$1]="+"$1
	}
	if (! neigh[$4])
		neigh[$4]="+"$4
	
	if ($8==0)
	{
		neigh[$1]=neigh[$1]" +"$4
		neigh[$4]=neigh[$4]" +"$1
	}
	else
	{
		neigh[$1]=neigh[$1]" -"$4
		neigh[$4]=neigh[$4]" -"$1
	}
}
END{

	for (str in neigh)
	{
		merge(str)
	}
	
	for (str in neigh)
	{
		test=1	# ready to print by default
		n=split(neigh[str], tmp, "+||-|| ")	# check for conflicts in orientation
		delete c
		for(i=1;i<=n;i++)
		{
			if (tmp[i]!="")
				c[tmp[i]]+=1
			if(c[tmp[i]]>1)
			{
				print ":( Some pairwise alignments are in conflict. Skipping merge block "neigh[str] "!" > "/dev/stderr"
				test=0
				break
			}
		}
		if (test && neigh[str]!="")
			print neigh[str]
	}
}