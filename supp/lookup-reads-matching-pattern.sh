#!/bin/bash

#### Description: Script to look for mapping positions of Hi-C reads that match a sequence pattern.
#### Usage: bash ./lookup-reads-matching-pattern.sh <sequence_pattern> <path_to_mnd_file>.
#### Input: string representing a pattern of interest and merged_nodups file. For the mnd file long format is necessary.
#### Output: stdout in wig format.
#### Parameters: none [TODO: k, s, bin can be made into paramteres].
#### Dependencies: Java; GNU Parallel.
#### Written by: Olga Dudchenko, version date 05/26/2020

pattern=$1
k=20;
s=5;
bin=100000;

export pattern
export k
export s
export bin

myfunc() {

awk -v pattern=$pattern -v k=$k -v s=$s -v b=${bin} -f <(cat - <<-'EOD'
function rc(str,	tmp,x, i){
	for(i=length(str); i>0; i--)
	{
		switch(substr(str,i,1))
		{
			case "A": tmp="T"
			break
			case "T": tmp="A"
			break
			case "G": tmp="C"
			break
			case "C": tmp="G"
			break
			default:
			tmp="N"
		}
		x=x""tmp;
	}
	return x;	
}
BEGIN{
	if(!k){k=40}
	if(!s){s=20}
	if(!b){b=1000000}
	if(length(pattern)<k)
	{
		if(pattern!=rc(pattern))
			test=pattern"|"rc(pattern)
		else
			test=pattern
	}
	else
	{
		for(i=1; i<=length(pattern)-k+s; i+=s)
		{
			testarr[substr(pattern,i,k)]=1
			testarr[rc(substr(pattern,i,k))]=1
		}
		for (i in testarr)
			test=test"|"i
		test=substr(test,2)
	}
#	print test > "/dev/stderr"
}
$0~test{
	if ($9>mapq) c[$2" "int($3/b)]++
	if ($12>mapq) c[$6" "int($7/b)]++
}
# $11~test{
# 	if($12==0)
# 		next
# 	if ($14~test)
# 		next
# 	c[$6" "int($7/b)]++
# 	next
# }
# $14~test{
# 	if($9==0)
# 		next
# 	if($11~test)
# 		next
# 	c[$2" "int($3/b)]++
# 	next
# }
END{
	for(i in c)
	{
		split(i,a," ")
		print a[1],a[2]*b, c[i]
	}
}
EOD
) $1
}

export -f myfunc

parallel -a ${2} --pipepart -j 80% --will-cite --block 1G --env myfunc myfunc | awk '{c[$1" "$2]+=$3}END{for(i in c){split(i,a," "); print a[1], a[2], c[i]}}' | sort -k 1,1 -k 2,2n | awk -v bin=$bin '$1!=prev{print "variableStep chrom="$1" span="bin; prev=$1}{print $2, $3}'
