## Script to plot coverage distribution to help qc coverage-related issues
## INPUT: coverage_wide wig file 
## OUTPUT: text histogram stdout
## Options: width, for width to bin the coverage values (default is 0.05); width and height (default 50 and 100 accordingly); start_interval, end_interval for dumping distribution only for a specific interval rather than the whole map.
## Written by: OD on 180928
BEGIN{
	if(!bin){bin=.05}
	if(!width){width=100}	# by default will build from 0X to 5X in terms of average coverage
	if(!height){height=50}
}
NR==1{
	if($0!~/^fixedStep/){print "Unrecognized file format. Exiting!" > "/dev/stderr"; exit}
	split($0,a,"=")
	res=a[length(a)]
	next
}
start_interval||end_interval{
        counter++
        if(counter*res<start_interval){next}
	if(counter*res>end_interval){exit}
}
{
	c[int($1/bin)]++
}
END{
	if(!length(c)){exit}
	
	for(i in c){
		if(c[i]>max){max=c[i]}
	}
	n=asorti(c, sortedc, "@ind_num_asc")
	
	print "" # leave one row up for ethical view
	for(j=1;j<=height;j++){
		str="  " # leave one left for axis
		k=1
		for(i=1;i<=width;i++){
			if(sortedc[k]!=(i-1)){str=str" ";continue}
			if((c[sortedc[k]]/max*height)>height-j+1){str=str"*"}else{str=str" "}
			k++
		}
		print str
	}
	str="0X"
	for(i=1;i<=width;i++){str=str"*"}
	str=str""(width*bin)"X"
	print str
}
