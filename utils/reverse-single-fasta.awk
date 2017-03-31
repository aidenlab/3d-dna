#### Description: Script to reverse complement a sequence.
#### USAGE: awk -f reverse-single-fasta.awk <path-to-single-fasta-file>
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
}
$0~/>/{$1=">RC_"substr($1,2); print $1; next}
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
