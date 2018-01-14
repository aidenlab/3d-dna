#### Description: Script to apply edits encoded in a 2D annotation file.
#### Usage: awk -v label1=<fraglabel> -v label2=<annolabel> -f edit-cprops-for-misassembled-contigs.awk <input-annotation-file> <cprops-file>
#### Input: Cprops file, 2D annotation file (0-based).
#### Output: New cprops in standard 3 column space separated format.
#### TODO (maybe): second label can be grabbed from the annotation file
#### Version: latest, after upstream modifications to annotation coordinate system.
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu, version 08/29/2016
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
# read in the annotation file: misassembled regions in original contigs only
{
	if (FILENAME==ARGV[1])
	{
		if (FNR!=1)
		{
			$2 += 1 # switch from 0-based to 1-based coordinate system in annotations
			if ($2 != 1)	# add preceeding nonmisassembled contig
			{
				if (blockcount[$1] == 0)
				{
					blockcount[$1]+=1
					start[$1" "blockcount[$1]]=1
					end[$1" "blockcount[$1]]=$2-1
				}
				else if (end[$1" "blockcount[$1]] != $2-1)
				{
					blockcount[$1]+=1
					start[$1" "blockcount[$1]]=end[$1" "(blockcount[$1]-1)]+1
					end[$1" "blockcount[$1]]=$2-1				
				}
			}
			blockcount[$1]+=1
			start[$1" "blockcount[$1]] = $2
			end[$1" "blockcount[$1]] = $3
			misasm[$1" "blockcount[$1]] = 1
		}
		next
	}
}
# read in the cprops file
{
	if (!($1 in blockcount))
	{
		counter+=1
		print $1, counter, $3
	}
	else
	{		
		for (i = 1; i <= blockcount[$1]; i++)
		{
			newname=$1""label1""i
			if (misasm[$1" "i])
			{
				newname=newname""label2
				if (end[$1" "i]-start[$1" "i]+1 == $3)	# whole scaffold annotation
					newname=$1""label2
			}
			counter+=1
			print newname, counter, end[$1" "i]-start[$1" "i] + 1
		}
		if (end[$1" "blockcount[$1]] != $3)
		{
			newname=$1""label1""(blockcount[$1]+1)
			counter+=1
			print newname, counter, $3-end[$1" "blockcount[$1]]
		}
	}
}