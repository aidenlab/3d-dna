#### Description: Awk script to calculate pth percentile.
#### Usage: awk -v p=<percentile> -f compute-centile.awk <sorted_dataset>
#### Input: Sorted dataset (num sort smaller to larger entry).
#### Parameters: p=<percentile>.
#### Output: stdout of pth centile.
#### Note: Implements linear interpolation between closest ranks method
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu. Version date 12/02/2016.
BEGIN{if(complement){p=100-p}}
{a[NR]=$1}
END{
	if (NR==1)
	{
		print a[1]
		exit
	}
	pos=p*(NR+1)/100
	if (pos<1)
		print a[1]
	else if (pos>=NR)
		print a[NR]
	else
	{
		d=pos-int(pos)
		print a[int(pos)]+d*(a[int(pos)+1]-a[int(pos)])
	}
}
