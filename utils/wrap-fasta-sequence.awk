#### Description: Script to wrap randomly formatter fasta files.
#### USAGE: awk -f wrap-fasta-sequence.awk <path-to-fasta-file>
#### Output: wrapped fasta-formatted stdout
#### Notes: Default is 80.

BEGIN{
	if (! block)
	{
		block=80
	}
}
{
	if ($0~/>/)
	{
		if (str != "")
			print str
		str=""
		counter=0
		i=0
		print
		next
	}
	str=str""$0
        counter+=length
	while (counter>=block)
	{
		print substr(str, i+1, block)
#		str=substr(str, block+1)
		counter-=block
		i+=block
	}
	str=substr(str, i+1)
	i=0
}
END{
        if (str != "")
        	print str
}
