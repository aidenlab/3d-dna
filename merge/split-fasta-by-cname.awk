## Script to split fasta by name to easily address individual contigs/scaffolds. Probably will be substituted by indexing later on.
## USAGE: awk -f split-fasta-by-name.awk <path_to_cprops> <path_to_fasta>
## Output: massive amount of individual fasta files.
## NOTE: Make sure the expected number of files does not exceed that allowed by the OS in a single folder...
FILENAME==ARGV[1]{
	cname[$1]=$2
}
$0~/>/{
	if(FNR!=1)
		close(filename)	
	filename=cname[substr($1,2)]
	print $1 > filename".fa"
	next
}
{
	print $0 > filename".fa"
}
