## Helper script to convert existing vcfs to psf and assembly files.
## USAGE: awk [ -v sample_name="NA12878" -v output_prefix="NA12878" -v chr=chrname(s) -v exclude_chr=list_of_chr_to_exclude -v ignore_filter=0 -v verbose=0 -v ignore_psf=0 ] -f vcf-to-assembly.awk <path_to_vcf_file> 
## INPUT: vcf file
## PARAMETERS: chr: can be a chromosome or a list of chromosomes, e.g. 1, chr1, 1|2, chr1|chr3|chr5 etc. TODO: more testing.
## Written by O.D. (olga.dudchenko@bcm.edu)
## Version: 5/21/19
## Modified: 11/19/20
## NOTE: requires testing with more crazy vcfs.
## TODO: deal with complex variants, maybe report how many SNPs are taken into consideration etc.

BEGIN{
	#defaults
	FS="\t"
	if(!length(sample_name)){sample_field=10}	# in case multisample vcf
	if(!length(output_prefix)){"basename "ARGV[1]" .vcf" | getline output_prefix}
# 	if(!length(output_prefix)){output_prefix="in"}
# 	if(length(chr)){output_prefix=output_prefix"."chr}
	if(length(chr)){chr="^("chr")$"}
	if(length(exclude_chr)){exclude_chr="^("exclude_chr")$"}
}

# skip metainfo lines
$0~/^##/{next}

# figure out which samples to parse out, if necessary
$0~/^#/{
	if (length(sample_name))
	{
		i=10 
		while( $i!=sample_name && i<=NF ){i++};
		if(i==NF+1){print ":| Warning! No "sample_name" sample found in VCF file. Exiting!" > "/dev/stderr"; force_exit=1; exit}
		sample_field=i;
	}; next
}
# if chr argument is passed ignore rows with different CHROM

(length(chr) && $1!~chr){next}

# if exclude_chr argument is passed ignore rows matching excluded CHROM
(length(exclude_chr) && $1~exclude_chr){next}

# check that there is expected number of columns
NF<sample_field{if(verbose){print ":| Unexpected number of fields, not sure how to parse VCF file. Ignoring VCF entry:" > "/dev/stderr"; print $0 > "/dev/stderr"}; next}

# check that REF is a simple variant
($4!~/^[ATCG]$/){if(verbose){print ":| Warning! Complex REF variant. Skipping VCF entry:" > "/dev/stderr"; print $0 > "/dev/stderr"}; next}

# check that ALT is a simple variant
($5!~/^[ATCG]$/)&&($5!~/^[ATCG],[ATCG]$/){if(verbose){print ":| Warning! Complex ALT variant. Skipping VCF entry:" > "/dev/stderr"; print $0 > "/dev/stderr"}; next}

# check that variant passes filter
(!ignore_filter && $7!="PASS"){if(verbose){print ":| Warning! FILTER entry is not PASS. Ignoring VCF entry:" > "/dev/stderr"; print $0 > "/dev/stderr"}; next}

# parse format and values for GT and PS
{	
	# FORMAT
	n=split($9, f, ":")
	# just in case
	if(f[1]!="GT"){if(verbose){print ":| Warning! Unexpected entries in FORMAT field. Ignoring VCF entry:" > "/dev/stderr"; print $0 > "/dev/stderr"}; next}

	# VALUES
	m=split($sample_field, a, ":")
	# just in case
	if(n!=m){
		if(verbose){print ":| Warning! Mismatch between FORMAT and values. Ignoring VCF entry:" > "/dev/stderr"; print $0 > "/dev/stderr"}; next}
	
	# GENOTYPE
	if(split(a[1],z,"\\|||/")!=2)
	{
		if(verbose){print ":| Warning! Unexpected genotype value. Ignoring VCF entry:" > "/dev/stderr"; print $0 > "/dev/stderr"}; next
	}
	if(z[1]==z[2])
	{
		next
	}
	
	split($5, var, ",")
	var[0]=$4
	
	if (!var[z[1]] || !var[z[2]]){if(verbose){print ":| Warning! GT inconsistent with variants listed. Ignoring VCF entry:" > "/dev/stderr"; print $0 > "/dev/stderr"}; next}
	
	var1=$1":"$2":"var[z[1]]
	var2=$1":"$2":"var[z[2]]
	
	counter+=2
	cprops[counter-1]=">"var1" "(counter-1)" "1000
	cprops[counter]=">"var2" "counter" "1000
	
	blockcprops[counter/2]=">"$1" "$2" "var[z[1]]" "var[z[2]]" "(counter/2)
	
	if (a[1]~/\//)
	{
		# singleton variant
		asm[counter-1]=" "counter-1
		asm[counter]=" "counter
		
		blockasm[counter/2]=" "(counter/2)
	}
	else
	{
		# prephased variant
		prephaseid=$1":-1"	# no PS field, interpret as belonging to the same set
		
		for(i in f)
		{
			if (f[i]=="PS"){prephaseid=$1":"a[i]}
		}
		if(!(prephaseid in block)){block[prephaseid]=counter}
		
		asm[block[prephaseid]-1]=asm[block[prephaseid]-1]" "(counter-1)
		asm[block[prephaseid]]=asm[block[prephaseid]]" "counter
		
		blockasm[block[prephaseid]/2]=blockasm[block[prephaseid]/2]" "(counter/2)
	}
}
END{
	# enforce exit
	if(force_exit){exit}

	for(i=1;i<=counter; i++)
	{
		print cprops[i] >output_prefix".assembly"
	}
	for(i=1;i<=counter; i++)
	{
		if(asm[i])
			print substr(asm[i],2) >output_prefix".assembly"
	}
	
	if (!ignore_psf)
		{
		for(i=1;i<=counter/2; i++)
		{
			print blockcprops[i] >output_prefix".psf" 
		}
		for(i=1;i<=counter/2; i++)
		{
			if (blockasm[i])
				print substr(blockasm[i],2) >output_prefix".psf" 
		}
	}
}
