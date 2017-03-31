#### Description: Script to apply edits encoded in a 2D annotation file to the mnd file.
#### Usage: awk -f edit-mnd-according-to-new-cprops.awk <edited-cprops-file> <old-mnd-file>
#### Input: Edited cprops file, mnd file.
#### Output: New mnd STDOUT in standard format.
#### WARNING: The cprops file is assumed to be consistent with the old mnd file. Very little internal checking is done!
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu on 08/29/2016; current version dated 12/10/16.
BEGIN{
	if (label1"")	# set consecutive piece counter label, e.g. fragment
	{
		label1=label1
	}
	else
	{
		print ":| Warning: No input for label1 was provided. Default for label1 is \":::fragment_\"" > "/dev/stderr"
		label1=":::fragment_"
	}
	if (label2"")	# set annotation piece label
	{
		label2=label2
	}
	else
	{
		print ":| Warning: No input for label2 was provided. Default for label2 is \":::debris\"" > "/dev/stderr"
		label2=":::debris"
	}
}
# read in the cprops file
{
	if (FILENAME==ARGV[1])
	{
		n=split($1, a , ""label1"||"label2"")
		if (n > 1)	# original scaffold was renamed
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
	# compact tentative output
	$10="-"; $11="-"; $13="-"; $14="-"; $15="-"; $16="-";

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
#		if ($3 < contigshift[$2" "i]+1 || $3 > shift[$2])
#		{
#			print "Annotations inconsistent with mnd. Exiting!"
#			print $0
#			exit 1
#		}

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
#		if ($7 < contigshift[$6" "i]+1 || $7 > shift[$6])
#		{
#			print "Annotations inconsistent with mnd. Exiting!"
#			print $0
#			exit 1
#		}

		$7=$7-contigshift[$6" "i]

		$6=contigname[$6" "i]
	}
	
	print	
}
