#### Description: Helper script to put unattempted back in after JB4A review if they were removed.
#### USAGE: awk -f unbundle-unattempted.awk <path-to-bundled-assembly> <path-to-unattempted-assembly>
#### Input: bundled and unattempted .assembly files
#### Output: .assembly-formatted stdout
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu on 171111.

# read in bundled assembly file
FILENAME==ARGV[1]{
	if ($0~/^>/){
		bundled_cprops[FNR]=$0
		cprops_counter=FNR
		next
	}
	bundled_asm[FNR-cprops_counter]=$0
	asm_counter=FNR-cprops_counter
	next
}
# read in unattempted assembly file
FILENAME==ARGV[2]{
	if ($0~/^>/){
		unattempted_cprops[FNR]=$0
		unattempted_cprops_counter=FNR
		next
	}
	if (NF!=1 || $1<0 || (FNR-unattempted_cprops_counter>1&&$1!=unattempted_asm[FNR-unattempted_cprops_counter-1]+1)){
		print "Unexpected format encountered in unattempted.asm. Exiting!" > "/dev/stderr"
		exit
	}
	unattempted_asm[FNR-unattempted_cprops_counter]=$0
}
END{
	split(bundled_cprops[cprops_counter], a, " ")
	if (a[1]!=">unattempted"){
		print "Unexpected format encountered in bundled.cprops. Exiting!" > "/dev/stderr"
		exit
	}	
	for (i in bundled_asm){
		if (bundled_asm[i]==a[2]){test=a[2]}
	}
	if (!test){
		print "Unexpected format encountered in bundled.asm. Exiting!" > "/dev/stderr"
		exit
	}
	
# print cprops
	for (i=1; i<=cprops_counter-1; i++) {
		print bundled_cprops[i]
	}
	for (i=1; i<=length(unattempted_cprops); i++){
		split(unattempted_cprops[i], a, " ")
		print a[1], cprops_counter, a[3]
		cprops_counter++
	}

# print asm
	for (i=1; i<=asm_counter; i++){
		if (bundled_asm[i]!=test){
			print bundled_asm[i] 
		}
	}
	shift=""
	for (i=1; i<=length(unattempted_asm); i++){
		if (shift==""){shift=test-unattempted_asm[i]}
		print unattempted_asm[i]+shift
	}
	
}
