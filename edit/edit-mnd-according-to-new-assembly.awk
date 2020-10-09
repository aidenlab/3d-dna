#### Description: Script to apply edits encoded in the assembly file to mnd file.
#### Usage: awk -f edit-mnd-according-to-new-assembly.awk <new-assembly-file> <old-mnd-file>
#### Input: New assembly file, old mnd file.
#### Output: New mnd STDOUT in standard format.
#### Parameters: label [string]; detailed [0/1]. 
#### WARNING: The assembly file is assumed to be consistent with the old mnd file. Very little internal checking is done!
#### WARNING: This script uses ":::" as default separator between draft sequence pieces.
#### WARNING: May loose a few contacts because of mnd shifting alignment position by S count in CIGAR.
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu on 08/29/2016; current version dated 12/10/16.
BEGIN{
	if (label"")	# sets consecutive piece counter label, by default :::
	{
		label=label
	}
	else
	{
		label=":::"
	}
}
# read in the cprops file
{
	if (FILENAME==ARGV[1])
	{
		if ($0!~/^>/){nextfile}
		$1=substr($1,2)
		n=split($1, a , ""label"")
		if (n > 1)	# original scaffold has changed
		{		
			counter[a[1]]+=1
			contigname[a[1]" "counter[a[1]]]=$1
			contigshift[a[1]" "counter[a[1]]]=shift[a[1]]
			shift[a[1]]+=$3
		}		
	next
	}
}
# mnd file: main loop
$0!=""{
	# compact tentative output unless explicitly requested not to
	if(!verbose)
	{
		$10="-"; $11="-"; $13="-"; $14="-"; $15="-"; $16="-";
	}

	if ((!($2 in counter)) && (!($6 in counter)))
	{
		print
		next
	}
	if ($2 in counter)
	{
		i=1	# identify the fragment to which to assign the contig
		while (contigshift[$2" "(i+1)] < $3 && i < counter[$2])
		{
			i++
		}

		$3 = $3-contigshift[$2" "i]
		
		$2=contigname[$2" "i]
	}	
	if ($6 in counter)
	{
		i=1
		while ( contigshift[$6" "(i+1)] < $7 && i < counter[$6])
		{
			i++
		}

		$7=$7-contigshift[$6" "i]

		$6=contigname[$6" "i]
	}
	
	print	
}
