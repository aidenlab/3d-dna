## Grab prefix or postfix of a given single fasta sequence
## Usage: awk -v test=<0/1> -v pos=<pos> -f grab-fix-sequence.awk <path-to-fasta>
## Parameters: test=0 if grabbing postfix; test=1 if grabbing prefix. pos represents the first (if postfix) or last (if prefix) base position. 
## Written by: OD
$0~/>/{

	if (test)
	{
        type="prefix"
	}
	if (!test)
	{
		type="postfix"
	}
	print $1":"type
	next
}
{
	c+=length
	if (c>=pos)
	{
		if (type=="postfix")
		{
			$0=substr($0,length-c+pos)
			test=1
		}
		else if (type=="prefix")
		{
			print substr($0,1,pos-c+length)
			exit
		}
	}
}
test
