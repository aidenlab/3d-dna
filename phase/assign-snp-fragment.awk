## Single-threaded script to assign snp fragments to reads for dangling SNP maps.
## USAGE: awk -f assign-snp-fragment.awk <path_to_snp_fragment_file> <path_to_mnd_file>
## SNP fragment file is a file formatted like restriction_site file, except it lists positions in between consecutive SNPs of interest. 
## Written by OD, version 210309

# read restriction_site_file
FILENAME==ARGV[1]{
	for(i=2;i<=NF;i++){
		position[$1" "(i-1)]=$i
	}
	counter[$1]=NF-1
	next
}
# read merged_nodups.txt file
{
	# do for read 1 and read 2 in read pair:
	for (half=1; half<=2; half++)
	{
		if (half==1)
		{
			chr=$2
			pos=$3
		}else{
			chr=$6
			pos=$7
		}
		
		upstream=1
		downstream=counter[chr]

	 	if (downstream==0 || pos <= position[chr" "upstream])
		{
			fragment[half]=0
			continue
		}

		if(pos > position[chr" "downstream])
		{
			fragment[half]=downstream
			continue
		}
	
		while (1)
		{
			tmp = int ((upstream+downstream)/2)
			if (pos <= position[chr" "tmp])
				downstream = tmp
			else if (pos > position[chr" "tmp])
				upstream = tmp

			if (upstream+1 == downstream)
			{
				break;
			}
		}
		fragment[half]=upstream
	}
	$4=fragment[1]
	$8=fragment[2]
	print
}
