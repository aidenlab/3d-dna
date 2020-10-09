## Generate an annotation track for merged sequences to help survey iterative merging results. Only called if unmerged annotations are passed to wrapper script.
## Color scheme adopted: white for secondary contigs; red for primary contigs wo upstream overlap; yellow for primary contigs with upstream overlap
## Written by: OD
BEGIN{
	FS="\t"
	OFS="\t"
	cc_intra="255,255,255"
	cc_break="255,0,0"
	cc_overlap="255,255,0"
}
# read in tiled annotation file
{
	if (FILENAME==ARGV[1])
	{
		if (FNR==1)
			next
		split($8,a," ")
		contigname=a[2]
		contigname=substr(contigname,3,length(contigname)-3)
		
		unmerged_left[contigname]=$2
		next
	}
}
# read in merger script annotation file
{
	if (FNR==1)
	{
		print
		next
	}
	contigname=$8
	if($7==cc_break)
	{
		last_left_merged = $2
		last_left_unmerged = unmerged_left[contigname]
		print "assembly", last_left_unmerged, last_left_unmerged+$3-$2, "assembly", last_left_unmerged, last_left_unmerged+$3-$2, cc_break, $8, last_left_unmerged, last_left_unmerged+$3-$2, last_left_unmerged, last_left_unmerged+$3-$2
	}
	if ($7==cc_intra)
	{
		print "assembly", last_left_unmerged+$2-last_left_merged, last_left_unmerged+$3-last_left_merged, "assembly",  unmerged_left[contigname], unmerged_left[contigname]+$3-$2, cc_intra, $8, last_left_unmerged+$2-last_left_merged, last_left_unmerged+$3-last_left_merged,  unmerged_left[contigname], unmerged_left[contigname]+$3-$2
	}
	if ($7==cc_overlap)
	{

		print "assembly", last_left_unmerged+$2-last_left_merged, last_left_unmerged+$3-last_left_merged, "assembly",  unmerged_left[contigname], unmerged_left[contigname]+$3-$2, cc_overlap, $8, last_left_unmerged+$2-last_left_merged, last_left_unmerged+$3-last_left_merged,  unmerged_left[contigname], unmerged_left[contigname]+$3-$2
		last_left_merged = $2
		last_left_unmerged = unmerged_left[contigname]
	}
}