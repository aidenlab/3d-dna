#### Description: Utility script to generate .assembly file from fasta. The order of scaffolds in the .assembly file matches that in the fasta file.
#### USAGE: awk -f generate-assembly-file-from-fasta.awk <path-to-fasta>
#### Output: assembly-formatted stdout

{
if ($0 ~ />/)
{
        if ( FNR != 1 )
		print ">"old_id, new_id, clength;
	old_id = substr($1,2)
	new_id += 1
	clength = 0
}
else
         clength += length
}
END{
	print ">"old_id, new_id, clength;
	for(i=1;i<=new_id; i++){
		print i
	}
}
