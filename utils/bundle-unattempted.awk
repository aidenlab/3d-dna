#### Description: Alternative helper script to bundle unattempted sequences at the end of the assembly. Currently used only when lifting maf files.
#### USAGE: awk [ -v input_size=${input_size} ] -f bundle-unattempted.awk <path-to-assembly>
#### Input: .assembly file
#### Output: two .assembly files with "bundled" and "unattempted" suffixes. Don't forget to unbundle!
#### NOTE: if no input_size is provided will try to bundle as much as possible?
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu on 200529.

BEGIN{start=0}
$0~/^>/{storecprops[$2]=$0; len+=$3; if(input_size && $3<input_size){input_start=$2; input_size=0}next}
{counter++; storeasm[counter]=$0}
NF==1&&$1>0{
	if(start && $1==prev+1)
		prev=$1
	else
	{
		start=$1; prev=$1
	}
	next
}
{start=0}
END{

	if(input_start)
	{
		if(input_start>=start){start=input_start}
		else
			print ":| Warning! Can't bundle everything!" > "/dev/stderr"
	}

	n=split(FILENAME,a,"/"); filename=substr(a[n],1,length(a[n])-8)
	
	bundledfilename=filename"bundled.assembly"
	unattemptedfilename=filename"unattempted.assembly"

	if(!start||start==$0){print ":| Warning: Nothing to bundle!" > "/dev/stderr"; system("ln -sf "FILENAME" "bundledfilename"; touch "unattemptedfilename); exit}
	
	if (start==1){print ":| Warning: Everything is in a single bundle!" > "/dev/stderr"}
	
	for(i=1;i<start;i++)
	{
		print storecprops[i]> bundledfilename
		split(storecprops[i],a," ")
		bundledlen+=a[3]
	}
	
	print ">unattempted", start, len-bundledlen > bundledfilename
	
	i=1; while(1)
	{
		if(storeasm[i]==start){break;}
		print storeasm[i]> bundledfilename
		i++
	}
	print start > bundledfilename
	
	for(i=start; i<=$0; i++)
	{
		print storecprops[i]> unattemptedfilename
	}
	for(i=start; i<=$0; i++)
	{
		print i> unattemptedfilename
	}
	
}