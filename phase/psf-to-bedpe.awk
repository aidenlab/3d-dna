## USAGE: awk -f psf-to-bedpe.awk <path_to_psf_file> 
## INPUT: psf file
## OUTPUT: stdout in bedpe format
## Written by O.D. (olga.dudchenko@bcm.edu)
## Version: 02/09/22

# print header
BEGIN{
    OFS="\t"
    lenPerSnp=1000
    end=0
    print "#chr1", "x1", "x2", "chr2", "y1", "y2", "name", "score", "strand1", "strand2", "color"
}
# read psf properties part
$0~/^>/{
	chr[$5]=substr($1,2)
	pos[$5]=$2
	var1[$5]=$3
	var2[$5]=$4
	next
}
!scale{
    scale=$(( 1 + 2*(NR-1)*lenPerSnp / 2100000000 ))
}
# read psf set part
{
    start=end; end+=int (NF*lenPerSnp/scale)
    print "assembly", start, end, "assembly", start, end, chr[$1]":"pos[$1]":"var1[$1], ".", ".", ".", "0,0,255" 
    start=end; end+=int (NF*lenPerSnp/scale)
    print "assembly", start, end, "assembly", start, end, chr[$1]":"pos[$1]":"var2[$1], ".", ".", ".", "0,0,255" 
}