#This file takes as input, the contig cprops, the asm file and a list of chromosome boundaries. We then output the first scaffold of each chromosome by locating the scaffold inside which this position lies.
BEGIN{
	if(ARGC!=5){print "USAGE: awk -f find_scaffolds_at_rabl_peak.awk <cprops> <asm> <chr_boundary_list> <scale_of_hic_map>"; exit -1;}
	global_pos = 1;
	scale = ARGV[4]
	
	ARGC=4;
}
function abs(x)
{
	if(x<0) return -x
	else return x
}
function sign(x)
{
	if(x<0) return -1
	else return 1
}
#Read in cprops
{
if(FILENAME==ARGV[1])
{
	cname[$1] = $2
	old_name[$2] = $1
	clen[$2] = $3
	next
}
}
#Read in asm
{
if(FILENAME==ARGV[2])
{
	n=split($0,a," ");
	for(i=1;i<=n;i++)
	{
		gpos[abs($i)] = global_pos
		scaf[global_pos] = abs($i)
		global_pos+=1

		ori[abs($i)] = sign($i)

		global_pos_bp[abs($i)] = pos_bp
		pos_bp += clen[abs($i)]
	}
	
	next
}
}
#Read in chromosome boundary list
{
#This file should have N_chr-1 lines, corresponding to the end positions of each chromosome, assuming that the hic map was built without dropouts
if(FILENAME==ARGV[3])
{
	if(FNR==1){First_Contig[FNR]=scaf[1];}

	chr_boundary[FNR] = $1*scale;

	#We locate the scaffold inside which this boundary lies
	for(i=1;i<=global_pos;i++)
	{
		if(global_pos_bp[scaf[i]]>chr_boundary[FNR]){break}
	}
	current_scaf = i-1;

	#print "Chr end for chromosome "FNR" lies at "chr_boundary[FNR]" inside scaffold "scaf[current_scaf] > "/dev/stderr"

	print scaf[current_scaf]


	First_Contig[FNR+1] = scaf[current_scaf]
	
	#print "First scaffold for chromosome "FNR+1" is: "First_Contig[FNR+1]#,max_index,current_max

	
	next
}
}

