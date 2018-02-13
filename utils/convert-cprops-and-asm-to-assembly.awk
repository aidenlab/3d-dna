#### Description: Helper script to convert from old cprops+asm format of storing assemblies to single-file .assembly format.
#### USAGE: awk -f convert-cprops-and-asm-to-assembly.awk <path-to-cprops> <path-to-asm>
#### Input: cprops and asm files
#### Output: xxx.assembly file based on xxx.asm file-name.
#### NOTE: will overwrite existing files!
#### NOTE: has some minimal user-proofness but could use more.
#### Written by: Olga Dudchenko- olga.dudchenko@bcm.edu. Version dated 180212.

# check that there is no obvious problem with the input
BEGIN{	
	if(substr(ARGV[1], length(ARGV[1])-6, length(ARGV[1]))!=".cprops" || substr(ARGV[2], length(ARGV[2])-3, length(ARGV[2]))!=".asm"){
		print ":( Unexpected file extensions. Check your input!" > "/dev/stderr"
		exit
	}
	n=split(substr(ARGV[1], 1, length(ARGV[1])-7), a, "/")
	m=split(substr(ARGV[2], 1, length(ARGV[2])-4), b, "/")
	if (a[n]!=b[m]){
		print ":| Warning, the input files have non-matching names. Make sure your inputs are internally compatible." > "/dev/stderr"
	}
}
# read in cprops
FILENAME==ARGV[1]{
	cprops_counter++
	print ">"$0 > substr(ARGV[2], 1, length(ARGV[2])-4)".assembly"
	next
}
{
	asm_counter+=NF
	print > substr(ARGV[2], 1, length(ARGV[2])-4)".assembly"
}
END{
	if (cprops_counter!=asm_counter){
		print ":( Input files not mutually compatible. Exiting!" > "/dev/stderr"
		system("rm "substr(ARGV[2], 1, length(ARGV[2])-4)".assembly")
	}
}