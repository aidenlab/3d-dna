#### Description: Script to preprocess reference data.
#### Input: Reference fasta file.
#### Output: Reference index stdout in format "sequence_name start_line_number end_line_number".
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu on 20190222
{
	if ($0~/^>/)
	{
		if ( NR != 1 )
			print seq_name, start_line, NR-1;
		seq_name = substr($1,2)
		start_line = NR
	}
}
END{
	print seq_name, start_line, NR
}