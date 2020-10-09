#### Description: Script to cut the fasta file to match the provided cprops file.
#### Usage: awk -v label1=<fraglabel> -v label2=<annolabel> -f edit-fasta-according-to-new-cprops.awk <cprops> <old-fasta-file>
#### WARNING: There are no internal checks for consistency. (Some will be added?)
#### Input: New cprops file, assumed to be compatible with the cprops file matching the 'old' fasta file and the fasta file itself. The cprops is assumed to be sorted with -k 1,1 -k 2,2n.
#### Output: Fasta STDOUT in standard format.
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu on 09/18/2016

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
	fasta_block_size=80 #wrap sequence
	i=1
}
# read in the new cprops file:
{
	if (FILENAME==ARGV[1])
	{
		split($1,a,""label1"||"label2"")
		editcounter[a[1]]+=1
		name[a[1]" "editcounter[a[1]]]=$1
		len[a[1]" "editcounter[a[1]]]=$3
		next
	}
}
# read in main fasta
{
	if ($0 ~ />/) # new contig or scaffold
	{
		# print leftovers from last contig or scaffold
		while (counter>=fasta_block_size)
		{
			print substr(str, 1, fasta_block_size)
			str=substr(str, fasta_block_size+1)
			counter-=fasta_block_size
		}
		if (str != "")
			print str
			
		if (current_subcontig_position != len[original_contig_name" "i])
		{
			system("echo \":( The length of fasta does not match that suggested by the cprops file. Exiting!\" >&2" )
			print 
			exit
		}
		
		original_contig_name=substr($1,2)
		i=1
		print ">"name[original_contig_name" "i]
		str=""
		counter=0	#string length counter for wrapping
		current_subcontig_position=0	#global string length since last >
		next
	}
		
	# wrap sequence	inside a block
	while (counter>=fasta_block_size)
	{
		print substr(str, 1, fasta_block_size)
		str=substr(str, fasta_block_size+1)
		counter-=fasta_block_size
	}
	
	# tentatively add next string
	str=str""$0
	counter+=length
	current_subcontig_position+=length

	if (current_subcontig_position <= len[contigname" "i])
	{
		next	#only wrap sequence
	}
	
	while (current_subcontig_position>len[original_contig_name" "i] && i < editcounter[original_contig_name]) #check when cross the end of current block
	{
#		print contigname, blockcount[contigname], start[contigname" "i], end[contigname" "i]
#		print length(str)-current_subcontig_position+end[contigname" "i]-1

		pre_str=substr(str,1,length(str)-current_subcontig_position+len[original_contig_name" "i])
		post_str=substr(str, length(str)-current_subcontig_position+len[original_contig_name" "i]+1)

		# finish printing sequence before the breakpoint
		str=pre_str
		counter-=length(post_str)
		while (counter>=fasta_block_size)	
		{
			print substr(str, 1, fasta_block_size)
			str=substr(str, fasta_block_size+1)
			counter-=fasta_block_size
		}
		if (str!="")
			print str
		
		# start writing next section
		i += 1
		print ">"name[original_contig_name" "i]
		str=post_str
		counter=length(str)
		current_subcontig_position=length(str)
	}
	
}
END{
	while (counter>=fasta_block_size)	# print last block
	{
		print substr(str, 1, fasta_block_size)
		str=substr(str, fasta_block_size+1)
		counter-=fasta_block_size
	}
	if (str!="")
		print str
}
