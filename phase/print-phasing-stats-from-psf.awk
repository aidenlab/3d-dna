## Helper script to dump everything from a psf file except for the largest phased block.
## USAGE: awk [ -v chr=chrname(s) ] -f extract-longest-blocks-from-psf.awk <path_to_psf_file>
## INPUT: psf file
## PARAMETERS: chr: can be a chromosome or a list of chromosomes, e.g. 1, chr1, 1|2, chr1|chr3|chr5 etc. TODO: more testing.
## Written by O.D. (olga.dudchenko@bcm.edu)
## Version: 1/27/21

BEGIN{
    if(length(chr)>0){
        split(chr,chromlist,"|")
        for(i in chromlist){
            c[chromlist[i]]=1
        }
    }
}
# create an output filename
# !outfile{ 
#     n = split(FILENAME, a, "/")
#     outfile=substr(a[n],1,length(a[n])-4)".largest.psf"
# }
# read properties bit
$0~/^>/{
    if((!chr)||(substr($1,2) in c))
    {
            chr_of_origin[$NF]=$1
            if(!var_counter[$1]){chrcounter++; chrname[chrcounter]=$1}
            var_counter[$1]++
    }    
}
# read phasing bit
!(chr_of_origin[$1] in var_counter){next}
NF>maxlength[chr_of_origin[$1]]{
    maxlength[chr_of_origin[$1]]=NF
    longest[chr_of_origin[$1]]=$0
}
END{
    for(i=1;i<=chrcounter;i++){
        #print longest[chrname[i]]>outfile
        printf "Phased on chromosome %s: %3.2f%% (%'d of %'d)\n", substr(chrname[i],2), maxlength[chrname[i]] / var_counter[chrname[i]] * 100, maxlength[chrname[i]], var_counter[chrname[i]]
    }
}