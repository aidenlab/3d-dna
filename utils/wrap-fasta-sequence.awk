#### Description: Script to wrap randomly formatter fasta files.
#### USAGE: awk -f wrap-fasta-sequence.awk <path-to-fasta-file>
#### Output: wrapped fasta-formatted stdout
#### Notes: Hardcoded 80, can make it a parameter
BEGIN{
	fasta_block_size=80
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
	while (counter>=fasta_block_size)
	{
		print substr(str, i+1, fasta_block_size)
#		str=substr(str, fasta_block_size+1)
		counter-=fasta_block_size
		i+=fasta_block_size
	}
	str=substr(str, i+1)
	i=0
}
END{
        if (str != "")
        	print str
}
