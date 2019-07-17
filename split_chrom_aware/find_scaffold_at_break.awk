#This file takes as input, the contig merged_nodups, the asm file and the cprops file for the draft scaffolds, and a list of chromosome boundaries. We then output the first scaffold of each chromosome by searching for a directionality index-esque score, but at the granularity of scaffolds instead of bp.
BEGIN{
	if(ARGC!=6){print "USAGE: awk -f find_scaffolds_at_break.awk <cprops> <asm> <merged_nodups> <chr_boundary_list> <scale_of_hic_map>"; exit -1;}
	global_pos = 1;
	pos_bp = 1;

	eps = 3.000001;

	radius = 4000000;#20Mb for directionality index

	cushion = 5000000;#5Mb for candidate search

	INT_MAX = 2140000000

	scale = ARGV[5]
	
	ARGC=5;
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
#Read contig merged_nodups and pre-process data
{
if(FILENAME==ARGV[3])
{
	c1 = cname[$2];
	c2 = cname[$6];
		
	#If both scaffolds are not in the asm then skip this line
	if(!( (c1 in gpos) && (c2 in gpos) )){next}


	#Consider only those read pairs with high mapq(>30)
	if(!( ($9>=30) && ($12>=30) )){next}

	if(ori[c1]==1)
	{
		offset_1 = $3;
	}
	else
	{
		offset_1 = clen[c1]-$3;
	}

	if(ori[c2]==1)
	{
		offset_2 = $7;
	}
	else
	{
		offset_2 = clen[c2]-$7;
	}
	
	#populate score entry for $2,$6, based on a radius (of 5Mb? Change?)
	n1 = global_pos_bp[c2]-global_pos_bp[c1]#+offset_2-offset_1

	if( (n1>0) && (n1<radius)){right_score[c1]+=1;left_score[c2]+=1;}
	else if( (n1<0) && (n1 > -1*radius) ){left_score[c1]+=1;right_score[c2]+=1;}

	next
}
}
#Read in chromosome boundary list
{
#This file should have N_chr-1 lines, corresponding to the end positions of each chromosome, assuming that the hic map was built without dropouts
if(FILENAME==ARGV[4])
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

	#Now we go cushion distance to the left and right of this scaffold and find the scaffold with the highest "ratio"(right_scort/left_score) 
	j=current_scaf;
	distance = 0;
	current_max = 0;
	max_index = 0; 
	while( (j>=1) && (distance <=cushion) )
	{
		ratio = (1.0*right_score[scaf[j]]+eps)/(1.0*left_score[scaf[j]]+eps);


		#if(scaf[j]=="21847"){print "j = "j" scaf = "scaf[j]" ratio = "ratio" current_max = "current_max" left_score = "left_score[scaf[j]]" right_score = "right_score[scaf[j]]}

                 if(ratio > current_max)
                 {
                         max_index = j;
                         current_max = ratio;
                 }
		#print "ratio = "ratio" max_index = "max_index" current_max = "current_max" j = "j" distance = "distance > "/dev/stderr"
		j = j-1;
		distance += clen[scaf[j]];
	}

	j = current_scaf+1;
	distance = 0;
	while( (j<global_pos) && (distance < cushion) )
	{
		ratio = (1.0*right_score[scaf[j]]+eps)/(1.0*left_score[scaf[j]]+eps);

		#if(scaf[j]=="21847"){print "j = "j" scaf = "scaf[j]" ratio = "ratio" current_max = "current_max" left_score = "left_score[scaf[j]]" right_score = "right_score[scaf[j]]}

		#if(scaf[j]=="26146"){print "j = "j" scaf = "scaf[j]" ratio = "ratio" current_max = "current_max" left_score = "left_score[scaf[j]]" right_score = "right_score[scaf[j]]}
		

		#if(scaf[j]=="14379"){print "j = "j" scaf = "scaf[j]" ratio = "ratio" current_max = "current_max" left_score = "left_score[scaf[j]]" right_score = "right_score[scaf[j]]}

                 if(ratio > current_max)
                 {
                         max_index = j;
                         current_max = ratio;
                 }


		#print "ratio = "ratio" max_index = "max_index" current_max = "current_max" j = "j" distance = "distance > "/dev/stderr"

		j = j+1;
		distance += clen[scaf[j]];
	}

	First_Contig[FNR+1] = scaf[max_index]
	
	print "First scaffold for chromosome "FNR+1" is: "First_Contig[FNR+1]#,max_index,current_max

	#print "Computed first scaffold for chromosome "FNR+1 > "/dev/stderr";
	
	next
}
}

