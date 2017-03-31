#This script accepts the original asm file(with the first line as the primary LIGer scaffold
#It also takes a list of the first scaffolds in each chromosome(starting from the second chromosome)
#It then splits the primary LIGer scaffold and outputs N_Chr lines
#Written by: Sanjit Batra
function abs(x)
{
	if(x<0){return -x}
	else{return x}
}
BEGIN{
n_chr=0;
current = 1;
}
#Read in the scaffold boundary names
{
if(FILENAME==ARGV[1])
{
	boundary[FNR]=$1
	n_chr+=1;
next
}
}
#Read in the asm file
{
if(FILENAME==ARGV[2])
{
	s = ""
	n=split($0,a," ");
	for(i=1;i<=n;i++)
	{
		if(abs($i)==boundary[current]){print s;s="";current+=1}
		s = s" "a[i]
		if(current>n_chr)
		{
			for(j=i+1;j<=n;j++){s=s" "a[j]}
			print s;
			exit +1
		}
	}
next
}
}
