#### Description: Awk script to calculate the depletion score along the genome. Can be used on either on the whole file or on a fragment of the dump file (e.g. together with GNU parallel).
#### Usage: awk -v bin_size=<bin_size> -v dep_size=<dep_size> -v sat_level=<sat_level> -f calculate-depletion-score.awk <hic_dump_file>
#### Input: Juicebox dump at <bin_size> resolution.
#### Parameters: bin_size: Juicebox dump resolution, dep_size: parameter that defines the size of the region across which the depletion score is averaged.
#### Output: stdout of depletion score in "pos score" format compatible with wig format. Note that the output is not sorted, except for the last position which is always listed last for convenience of parallelization. 
#### Note: In this version of the script the score for bin i is calculated across a triangular area between lines x=y=i and y-x=dep_size.
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu. Version date 11/30/2016.
$3!="NaN"{
		
	if ($2-$1>dep_size)	# ignore data away outside the region of interest
	{
		next
	}
	if ($3>=sat_level)	# saturate matrix
	{
		$3=sat_level
	}
	for (i=$1+bin_size; i<=$2-bin_size; i+=bin_size)	# contribute to scores
	{
		c[i]+=$3
	}
}
END{
	for(i in c)
	{
			print i, c[i]
	}
}
