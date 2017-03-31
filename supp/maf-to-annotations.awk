BEGIN{OFS="\t"; print "chr1", "x1", "x2", "chr2", "y1", "y2", "color", "id", "X1", "X2", "Y1", "Y2"}
#BEGIN{id1="drosophila_melanogaster"; id2="aedes_aegypti"}
#BEGIN{id1="AgamP4"; id2="AaegL3"}
$0!~/#/{

	n=split($2,tmp,".")
	id=tmp[1]
	tmpstr=tmp[2]
	for (i=3; i<=n; i++)
		tmpstr=tmpstr"."tmp[i]
	contig[id]=tmpstr
	orientation[id]=$5
	if (orientation[id]=="+")
		start_coordinate[id]=$3
	else
		start_coordinate[id]=$6-$3
	end_coordinate[id]=start_coordinate[id]+$4
	if (id==id2)
	{
		print contig[id2], start_coordinate[id2]-1, end_coordinate[id2]-1, contig[id2], start_coordinate[id2]-1, end_coordinate[id2]-1, "0,0,0", orientation[id2]""contig[id1]":"(start_coordinate[id1]-1)"-"(end_coordinate[id1]-1), start_coordinate[id2]-1, end_coordinate[id2]-1, start_coordinate[id2]-1, end_coordinate[id2]-1
	}
}
END{
		print contig[id2], start_coordinate[id2]-1, end_coordinate[id2]-1, contig[id2], start_coordinate[id2]-1, end_coordinate[id2]-1, "0,0,0", orientation[id2]""contig[id1]":"(start_coordinate[id1]-1)"-"(end_coordinate[id1]-1), start_coordinate[id2]-1, end_coordinate[id2]-1, start_coordinate[id2]-1, end_coordinate[id2]-1
}
