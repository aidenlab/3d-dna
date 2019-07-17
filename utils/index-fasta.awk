#### Description: Script to preprocess reference data. Designed to substitute two scripts: generate cprops and index fasta. Used for adding reference data to the map.
#### Input: Reference fasta file.
#### Output: Reference index stdout in format "chrom/scaf/contig_name start_byte chrom/scaf/contig_length".
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu on 06/14/2016
#### Notes: Could be parallelized
{
	c += length + 1

	if ($0~/>/)
	{
		if ( NR != 1 )
			print scaf_name, start_byte, clength;
		scaf_name = substr($1,2)
		start_byte = sprintf("%4i", c+1)
		clength = 0
	}
	else
		clength += length
}
END{
	print scaf_name, start_byte, clength
}
