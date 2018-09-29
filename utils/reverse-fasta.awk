#### Description: Script to reverse complement a sequence.
#### USAGE: awk -f reverse-fasta.awk <path-to-fasta-file>
#### Output: fasta-formatted stdout

function substitute(char){
## TODO: handle non-standard?
	if (char=="A")
		return "T"
	else if (char=="T")
		return "A"
	else if (char=="G")
		return "C"
	else if (char=="C")
		return "G"
	else if (char=="a")
		return "t"
	else if (char=="t")
		return "a"
	else if (char=="c")
		return "g"
	else if (char=="g")
		return "c"
	else if (char=="N")
		return "N"
	else if (char=="n")
		return "n"
	else
	{
		if(!warning)
		{
			print ":| WARNING: fasta contains ambiguous bases!" > "/dev/stderr"
			warning=1
		}
		if (char=="Y")
			return "R"
		else if (char=="y")
			return "r"
		else if (char=="R")
			return "Y"
		else if (char=="r")
			return "y"
		else if (char=="W")
			return "W"
		else if (char=="w")
			return "w"
		else if (char=="S")
			return "S"
		else if (char=="s")
			return "s"
		else if (char=="K")
			return "M"
		else if (char=="k")
			return "m"
		else if (char=="M")
			return "K"
		else if (char=="m")
			return "k"
		else if (char=="D")
			return "H"
		else if (char=="d")
			return "h"
		else if (char=="H")
			return "D"
		else if (char=="h")
			return "d"
		else if (char=="V")
			return "B"
		else if (char=="v")
			return "b"
		else if (char=="B")
			return "V"
		else if (char=="b")
			return "v"
		else if (char=="X")
			return "X"
		else if (char=="x")
			return "x"
		else {
			print "!ERROR: Unknown base, exiting!" > "/dev/stderr"
			exit
		}
	}
}
$0~/>/{
	for(i=n;i>=1;i--)
	{
		str=""
		k=split(line[i],a,"")
		for(s=k;s>=1;s--)
		{
			str=str""substitute(a[s])
		}
		print str
	}
	$1=">RC_"substr($1,2)
	print $1
	n=0
	next
}
{n++; line[n]=$0}
END{
for(i=n;i>=1;i--)
{
	str=""
	k=split(line[i],a,"")
	for(s=k;s>=1;s--)
	{
		str=str""substitute(a[s])
	}
	print str
}
}
