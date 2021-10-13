#### Description: Helper script to extract the cprops portion from an assembly assembly file.
#### USAGE: awk -f extract-cprops-from-assembly.awk <path-to-assembly>
#### Input: ".assembly" file.
#### Output: xxx.cprops file where xxx name is based on xxx.assembly file name.
#### NOTE: will overwrite existing files!
#### Written by: Olga Dudchenko- olga.dudchenko@bcm.edu. Version dated 210702.

$0~/^>/{
	$1=substr($1,2)
	print > substr(FILENAME, 1, length(FILENAME)-9)".cprops"
}