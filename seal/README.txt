change log:
—————————————————————
version september 22, 2018.
	+ fixed bug occurring with the last sequence in the properties file was among those sealed
—————————————————————
version april 21, 2018.
	+ rewrote to take and output .assembly format
	+ modified algorithm to correctly seal intermittent groups rather than just triples
	+ added handling of unattempted for easy bundling
	TODO: maybe will separate bundling size option from stringency options instead of using input size for both
—————————————————————
version may 18, 2017.
	+ rewrote to accomodate for single very large contigs/scaffolds. These are now not moved to the end of the scaffold list but retail their relative position in the asm file.
version oct 2017.	
	+ rewrote to account for a rare bug
—————————————————————
