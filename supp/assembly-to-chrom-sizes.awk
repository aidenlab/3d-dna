BEGIN{OFS="\t"}
$0~/^>/{
	len[$2]=$3;
	next
}
{counter++}
!n||counter<=n{
	gsub("-","")
	c=0
	for(i=1;i<=NF;i++){c+=len[$i]}
	print "HiC_scaffold_"counter, c
}
