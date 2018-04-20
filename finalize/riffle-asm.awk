## Helper script to riffle the assembly with a gap contig
## Written by: OD
{
str=$1
for(i=2;i<=NF;i++)
{
	str=str" "riffle" "$i
}
print str
}
