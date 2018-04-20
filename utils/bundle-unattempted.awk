#### Description: Helper script to bundle unattempted for JB4A review in case there are too many.
#### USAGE: awk -v input_size=${input_size} -f bundle-unattempted.awk <path-to-assembly>
#### Input: .assembly file
#### Output: two assembly files with "bundled" and "unattempted" suffixes. Use bundled for review, but don't forget to unbundle!
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu on 171111.

BEGIN{
	if(!input_size){
		input_size=15000 # defaults
	}
}
# read in cprops
$0~/^>/{
	cprops[FNR]=$0
	if($1!~/:::fragment_/ && $3<input_size){
		unattempted[$2]=1
	}
	counter++
	next
}
# read in asm
{
	scaffold[NR-counter]=$0
}
END{
	if (NF==1 && $1==counter && $1 in unattempted){
		last=$1
	}
	for(i=NR-counter-1; i>=1; i--){
		n=split(scaffold[i],a," ");
		if (n==1 && a[1]==last-1 && a[1] in unattempted){
			last=a[1]
		}else{
			break
		}
	}
	for (k=1; k<=last-1; k++){
		print cprops[k] > substr(FILENAME,1,length(FILENAME)-8)"bundled.assembly"
	}
	for (k=last; k<=counter; k++){
		split(cprops[k],a," ")
		pseudo+=a[3]
	}
	print ">unattempted", last, pseudo > substr(FILENAME,1,length(FILENAME)-8)"bundled.assembly"
	for (k=last; k<=counter; k++){
		print cprops[k] > substr(FILENAME,1,length(FILENAME)-8)"unattempted.assembly"
	}
	
	for (k=1; k<=i+1; k++){
		print scaffold[k] > substr(FILENAME,1,length(FILENAME)-8)"bundled.assembly"
	}
	for (k=i+1; k<=NR-counter; k++){
		print scaffold[k] > substr(FILENAME,1,length(FILENAME)-8)"unattempted.assembly"
	}
	
}