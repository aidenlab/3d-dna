## Construct tiling path through the scaffolds of each individual connected component/cluster
## Written by: OD
function calculate_leftmost(c1, c2, algnstr,	a){
	split(algnstr, a, "\t")
	if (c1>0 && c2>0)
	{
		if (oname[c1]==a[1])
			return pos[c1]+a[2]-a[5]
		else
			return pos[c1]+a[5]-a[2]
	}
	else if (c1>0 && c2<0)
	{
		if (oname[c1]==a[1])
		{
			return pos[c1]+a[2]-len[a[4]]+a[6]
		}
		else
		{
			return pos[c1]+a[5]-len[a[1]]+a[3]
		}
	}
	else if (c1<0 && c2>0)
	{
		if(oname[-c1]==a[1])
		{
			return pos[c1]+len[a[1]]-a[3]-a[5]
		}
		else
		{
			return pos[c1]+len[a[4]]-a[6]-a[2]
		}
	}
	else
	{
		if(oname[-c1]==a[1])
		{
			return pos[c1]+len[a[1]]-len[a[4]]-a[3]+a[6]
		}
		else
		{
			return pos[c1]+len[a[4]]-len[a[1]]-a[6]+a[3]
		}
	}
}
function DFS(element,	pair, b){
	for ( pair in align )
	{
		split(pair, b)
		if (b[1]==element || b[1]==-element)
		{
			if ((b[2] in pos) || (-b[2] in pos))
				continue
			tmp=align[b[1]" "b[2]" "b[3]]
#			delete align[b[1]" "b[2]" "b[3]]
#			delete align[b[2]" "b[1]" "b[3]]

			if (element>0&&b[3]==0||element<0&&b[3]==1)
			{
				pos[b[2]]=calculate_leftmost(element, b[2], tmp)
				DFS(b[2])
			}
			else
			{
				pos[-b[2]]=calculate_leftmost(element, -b[2], tmp)
				DFS(-b[2])
			}
		}
	}
}
## read in cprops
FILENAME==ARGV[1]{
	oname[$2]=$1
	oname[-$2]=$1
	cname[$1]=$2
	len[$1]=$3
	next
}
## read in alignments
FILENAME==ARGV[2]{
	align[cname[$1]" "cname[$4]" "$8]=$0
	align[cname[$4]" "cname[$1]" "$8]=$0
	next
}
# individual cluster asm
{
	delete pos
	pos[$1]=1
	DFS($1)

	n = asort(pos, dest)
	
	outstr=""
	for(i=1; i<=n; i++)
	{
		for(k in pos)
			if (pos[k]==dest[i])
			{
				outstr=outstr" "k
				delete pos[k]
			}
	}
	print substr(outstr,2)		
}