#This file takes as input the first row of the juicebox dump at 50kb resolution and a list of coarse chromosome boundary positions.
#It then searches for the maximum value in a radius of 1Mb centered at each of the above positions.
BEGIN{
	res = 50000;
}
#Read in Juicebox dump
{
if(FILENAME==ARGV[1])
{
	serial_num[$1] = FNR
	position_of_serial[FNR]=$1
	val[FNR] = $2
next
}	
}
#Read in positions
{
if(FILENAME==ARGV[2])
{
	pos = serial_num[$1]
	max = -1;
	for(i=pos-20;i<=pos+20;i++)
	{
		if(i in val)
		{
			if(max<val[i]){max = val[i];position[FNR]=position_of_serial[i+1]}
		}
	}
	num_lines+=1;
next
}	
}
END{
	for(j=1;j<=num_lines;j++)
	{
		print position[j];
	}
}
