#### Description: Helper script to convert from new assembly format to old cprops+asm format for storing assemblies.
#### USAGE: awk -f convert-assembly-to-cprops-and-asm.awk <path-to-assembly>
#### Input: ".assembly" file.
#### Output: xxx.cprops and xxx.asm files based on xxx.assembly file name.
#### NOTE: will overwrite existing files!
#### Written by: Olga Dudchenko- olga.dudchenko@bcm.edu. Version dated 180212.

$0~/^>/{
	$1=substr($1,2)
	print > substr(FILENAME, 1, length(FILENAME)-9)".cprops"
	next
}
{
	print > substr(FILENAME, 1, length(FILENAME)-9)".asm"
}