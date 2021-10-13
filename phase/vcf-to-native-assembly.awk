## Helper script to extract contig names and sizes from vcfs so that we don't need a chrom.sizes input or something.
## USAGE: awk [ -v chr=<chrname(s)> ] -f vcf-to-native-assembly.awk <path_to_vcf_file> 
## INPUT: vcf file
## PARAMETERS: chr: can be a chromosome or a list of chromosomes, e.g. 1, chr1, 1|2, chr1|chr3|chr5 etc.
## Written by O.D. (olga.dudchenko@bcm.edu)
## Version: 1/27/21
## NOTE: requires testing with more crazy vcfs.

BEGIN{
    FS="\t";
    if(length(chr)>0){
        split(chr,chromlist,"|")
        for(i in chromlist){
            c[chromlist[i]]=1
        }
    }
}
$0!~/^#/{exit}
$0~/^##contig=<ID=[0-z_.-]+,length=[0-9]+>$/{
    split(substr($0,14),a,",")
    name=a[1]
    split(a[2],b,"=")
    len=substr(b[2],1,length(b[2])-1)
    if(!chr||(name in c))
    {
        counter+=2
        print ">"name"-r", counter-1, len
        print ">"name"-a", counter, len
    }
}
END{
    if(!counter){print ":( No expected metadata found in the vcf file. Exiting!" > "/dev/stderr"; exit 1}
    if(counter>100){print ":| Warning: number of sequences to be included in chrom.sizes  exceeds 100." > "/dev/stderr"}
    for(i=1;i<=counter;i++){print i}
}
