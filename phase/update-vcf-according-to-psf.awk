## Helper script to update phasing information in a vcf file in accordance with a psf file.
## USAGE: awk -f update-vcf-using-psf.awk <path_to_psf_file> <path_to_vcf_file> 
## INPUT: psf file, vcf file
## Written by O.D. (olga.dudchenko@bcm.edu)
## Version: 3/2/21
## NOTE: With more testing, this script is supposed to substitute psf-to-vcf.awk
## WARNING: not tested on multisample vcf
## WARNING: It may make sense in the future to get rid of prephasing in the vcf file, e.g. associated with non-primitive variants etc. Right now these are kept with the original phasing block assignment, even if that block has been fused with something else. This is in line with other vcf processing tools that keep phasing block info the same even if other variants associated with the block are filtered out.
BEGIN{
	if(!length(sample_name)){sample_field=10}	# in case multisample vcf
    OFS="\t"
}
# read in psf
FILENAME==ARGV[1]{
    if($0~/^>/){
        id[substr($1,2)" "$2" "$3" "$4]=$5
        pos[$5]=$2
    }else if(NF!=1)
    {
    	ps=pos[$1]
        for(i=1;i<=NF;i++)
        {
            if($i>0)
                phase[$i]=ps
            else
                phase[-$i]=-ps			
        }
    }
    next
}

# read in vcf
{
	FS="\t"
}
# print metainfo lines
$0~/^##FORMAT=<ID=/{has_format_metainfo=1}
$0~/^##FORMAT=<ID=PS,/{has_ps_metainfo=1}
$0~/^##/{print; next}
has_format_metainfo&&(!has_ps_metainfo){print "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phasing ID information, where each unique ID connects records within a phasing group\">"; has_ps_metainfo=1}
$0~/^#/{
    print "##source=3D-DNA";
    print;
    # figure out which samples to parse out, if necessary
	if (length(sample_name))
	{
		i=10 
		while( $i!=sample_name && i<=NF ){i++};
		if(i==NF+1){print ":| Warning! No "sample_name" sample found in VCF file. Exiting!" > "/dev/stderr"; force_exit=1; exit}
		sample_field=i;
	}
    next
}

## THIS IS NOT NECESSARY AS THESE ARE NOT SUPPOSED TO BE IN THE PSF ANYWAY
# check that REF is a simple variant
($4!~/^[ATCG]$/){print; next}
# check that ALT is a simple variant
($5!~/^[ATCG]$/)&&($5!~/^[ATCG],[ATCG]$/){print; next}
# check that variant passes filter
(!ignore_filter && $7!="PASS"){print; next}

# parse format and values for GT and PS
{	
	# FORMAT
	n=split($9, f, ":")
	# just in case
	if(f[1]!="GT"){if(verbose){print ":| Warning! Unexpected entries in FORMAT field. Ignoring VCF entry:" > "/dev/stderr"; print $0 > "/dev/stderr"}; next}
    
    PSfield=0; for(i in f){if(f[i]=="PS"){PSfield=i}}

	# VALUES
	m=split($sample_field, a, ":")
	# just in case
	if(n!=m){if(verbose){print ":| Warning! Mismatch between FORMAT and values. Ignoring VCF entry:" > "/dev/stderr"; print $0 > "/dev/stderr"}; next}

	# GENOTYPE
	if(split(a[1],z,"\\|||/")!=2)
	{
		if(verbose){print ":| Warning! Unexpected genotype value. Ignoring VCF entry:" > "/dev/stderr"; print $0 > "/dev/stderr"}; next
	}

    # REMOVE OLD PHASING INFORMATION??? Sort of cannot confirm or deny?

	if(z[1]==z[2])
	{
		print; next
	}

    split($5, var, ",")
	var[0]=$4
	
	if (!var[z[1]] || !var[z[2]]){if(verbose){print ":| Warning! GT inconsistent with variants listed. Ignoring VCF entry:" > "/dev/stderr"; print $0 > "/dev/stderr"}; next}

    if(phase[id[$1" "$2" "var[z[2]]" "var[z[1]]]])
    {
        phase[id[$1" "$2" "var[z[1]]" "var[z[2]]]]=-phase[id[$1" "$2" "var[z[2]]" "var[z[1]]]]
    }
    
    if (phase[id[$1" "$2" "var[z[1]]" "var[z[2]]]])
    {
        if(phase[id[$1" "$2" "var[z[1]]" "var[z[2]]]]>0){
            a[1]=z[1]"|"z[2]
            gsub("^[^:]*:", a[1]":", $sample_field)
        }else{
            a[1]=z[2]"|"z[1]
            gsub("^[^:]*:", a[1]":", $sample_field)
            phase[id[$1" "$2" "var[z[1]]" "var[z[2]]]]=-phase[id[$1" "$2" "var[z[1]]" "var[z[2]]]]
        }
        
        if(PSfield){
            $sample_field=a[1]
            for(i=2;i<=m;i++){if(i==PSfield){a[i]=phase[id[$1" "$2" "var[z[1]]" "var[z[2]]]]}; $sample_field=$sample_field":"a[i]}
        }else{
            $9=$9":PS"
            $sample_field=$sample_field":"phase[id[$1" "$2" "var[z[1]]" "var[z[2]]]]
        }
    }

    print;

}