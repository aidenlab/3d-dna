BEGIN{
	OFS="\t";
	if (!target_file_name){
		target_file_name="target_annotations.txt"
	}
	if (!query_file_name){
		query_file_name="query_annotations.txt"
	}
	
	print "chr1", "x1", "x2", "chr2", "y1", "y2", "color", "id", "X1", "X2", "Y1", "Y2" > ""target_file_name""
	print "chr1", "x1", "x2", "chr2", "y1", "y2", "color", "id", "X1", "X2", "Y1", "Y2" > ""query_file_name""
}
!$0{next}
$0~/^a score=/{
	include = 0
	split($0,a,"=")
	score=a[2]
	if(score>score_threshold){
		include = 1
	}
	next
}
(include){
	contig[include]=$2
	if (id1 && include==1){split(contig[include],tmp,id1); contig[include]=tmp[2]}
	if (id2 && include==2){split(contig[include],tmp,id2); contig[include]=tmp[2]}

	orientation[include]=$5
	if (orientation[include]=="+")
		start_coordinate[include]=$3
	else
		start_coordinate[include]=$6-$3-$4
	
	end_coordinate[include]=start_coordinate[include]+$4
	
	if (include==2) {
		
		# target
		print contig[1], start_coordinate[1], end_coordinate[1], contig[1], start_coordinate[1], end_coordinate[1], "0,0,0", "+"contig[1]":"(start_coordinate[1])"-"(end_coordinate[1])":"score, start_coordinate[1], end_coordinate[1], start_coordinate[1], end_coordinate[1] > ""target_file_name""
		
		# query
		print contig[2], start_coordinate[2], end_coordinate[2], contig[2], start_coordinate[2], end_coordinate[2], "0,0,0", orientation[2]""contig[1]":"(start_coordinate[1])"-"(end_coordinate[1])":"score, start_coordinate[2], end_coordinate[2], start_coordinate[2], end_coordinate[2] > ""query_file_name""

		}
		
		include++	
	if (include > 3) {
 		print "Unrecognized format" > "/dev/stderr"
 	}
}
