#### Description: Script to map 2D annotation files from a given assembly to the original input scaffolds/contigs according to a particular assembly file.
#### Usage: awk -v scale=<scale> -v split_annotations=<0/1> -f lift-assembly-annotations-to-input-annotations.awk <assembly> <assembly-annotation-file>
#### Input: Assembly file and 2D annotation file. Note that annotations are expected to be in the bed-like 0-based coordinate format (0-1000; 1000-2000; 2000-3000), 12-column assembly tab separated format.
#### WARNING: Works for diagonal annotations only at this point. TODO!
#### Output: 2D annotation file STDOUT in assembly 12-column tab separated format.
#### Parameters: parameter to scale coordinates in final output, can be both >1 (to squeeze the map) or <1 (to stretch the map out, default scale = 1). Never needed, kept just in case for cases when individual scaffolds are > 2.1 Gb.
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu, version date 20180525

BEGIN{
	OFS="\t"
	if (!scale)
	{
		scale=1
	}
	pos = 1
}
function abs(value)
{
  return (value<0?-value:value)
}
#read in cprops
{
	if (FILENAME==ARGV[1] && $0~/^>/)
	{
		oname[$2] = substr($1,2)
		clength[$2] = $3
		next
	}
}
#read in asm
{
	if (FILENAME==ARGV[1])
	{
		n=split($0,a)

		for (i=1; i<=n; i++)
		{
			if (!( abs(a[i]) in clength ))
			{
				print ":( Assembly file does not match cprops file. Exiting!" > "/dev/stderr"
				exit 1
			}
			counter += 1
			asm[counter] = abs(a[i])
			if ( abs(a[i]) in globalpos )
			{
				print ":( Assembly file contains multiple entries of the same contig. Exiting!" > "/dev/stderr"
				exit 1
			}
			globalpos[abs(a[i])] = pos
			if (abs(a[i])!=a[i])
			{	
				orientation[abs(a[i])] = 16
			}
			else
			{
				orientation[abs(a[i])] = 1
			}	
			pos += clength[abs(a[i])]
		}
		next
	}
}
#read in annotation file
{FS="\t"}
FNR==1{print; next}
{
	remember=$0; splitstartx=$9; tryagainx=1; i=1;
	while(tryagainx)
	{

		$0=remember
	
		# axis 1: fields $9 and $10
	
		$9=splitstartx
		# lift unscaled coordinates, switch to 1-based coordinate system
		$2=$9+1
		$3=$10
	
		# figure out which contig. Maybe TODO: optimize by sorting the original annotation file

		while ( globalpos[asm[i]] + clength[asm[i]] -1 < $2 && i <= counter)
		{
			i++
		}
	
		if ( i > counter)
		{
			print ":( Annotations are our of range covered by suggested contigs. Exiting!" > "/dev/stderr"
			exit 1
		#next
		}

		# TODO: handle splitting of annotations here.	
		if ( $3 > globalpos[asm[i]] + clength[asm[i]] )
		{
			if(split_annotations)
			{
				$3 = globalpos[asm[i]] + clength[asm[i]] - 1
				splitstartx=$3
				tryagainx=1
			}else{
				print ":| Warning: Annotations span split contigs. Ignoring!" > "/dev/stderr"
				next
			}
		}
		else
		{
			tryagainx=0
		}
	
		$1 = asm[i]

		if ( orientation[$1] == 1)
		{
			tmp2 = $2 - globalpos[$1] + 1
			tmp3 = $3 - globalpos[$1] + 1
		
		}
		else
		{
			tmp2 = -$3 + globalpos[$1] + clength[$1]
			tmp3 = -$2 + globalpos[$1] + clength[$1]

	# 		if ($8~/^-/)
	# 		{
	# 			sub("^\\-","+", $8)
	# 			sub(" \\(\\-"," (+", $8)
	# 		}
	# 		else if ($8~/^+/)
	# 		{
	# 			sub("^\\+","-", $8)
	# 			sub(" \\(\\+"," (-", $8)
	# 		}
		}
		$1=oname[$1]
		$9=tmp2-1; $10=tmp3; 
		$2=int($9/scale); $3=int($10/scale);
		remembery=$0
		# set up for splitting y axis
		splitstarty=$11; tryagainy=1; j=1;	
		while (tryagainy){
			$0=remembery

			# axis 2: fields $11 and $12

			$11=splitstarty
			# lift unscaled coordinates, switch to 1-based coordinate system
			$5=$11+1
			$6=$12

			# figure out which contig. Maybe TODO: optimize by sorting the original annotation file

			while ( globalpos[asm[j]] + clength[asm[j]] -1 < $5 && j <= counter)
			{
				j++
			}

			if ( j > counter)
			{
				print ":( Annotations are our of range covered by suggested contigs. Exiting!" > "/dev/stderr"
				exit 1
			}

			# TODO: handle splitting of annotations here.	
			if ( $6 > globalpos[asm[j]] + clength[asm[j]] )
			{
				if(split_annotations)
				{
					$6 = globalpos[asm[j]] + clength[asm[j]] - 1
					splitstarty=$6
					tryagainy=1
				}else{
					print ":| Warning: Annotations span split contigs. Ignoring!" > "/dev/stderr"
					next
				}
			}
			else
			{
				tryagainy=0
			}

			$4 = asm[j]

			if ( orientation[$4] == 1)
			{
				tmp2 = $5 - globalpos[$4] + 1
				tmp3 = $6 - globalpos[$4] + 1
	
			}
			else
			{
				tmp2 = -$6 + globalpos[$4] + clength[$4]
				tmp3 = -$5 + globalpos[$4] + clength[$4]

	# 			if ($8~/^-/)
	# 			{
	# 				sub("^\\-","+", $8)
	# 				sub(" \\(\\-"," (+", $8)
	# 			}
	# 			else if ($8~/^+/)
	# 			{
	# 				sub("^\\+","-", $8)
	# 				sub(" \\(\\+"," (-", $8)
	# 			}
			}
			$4=oname[$4]
			$11=tmp2-1; $12=tmp3;
			$5=int($11/scale); $6=int($12/scale);
		
		
			# print output
			print		
		
		}
	}
}