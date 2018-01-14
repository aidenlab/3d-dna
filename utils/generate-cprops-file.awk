#### Description: Generate cprops file describing one-to-one correspondence between original draft sequences names and internal naming convention. Includes total lengths.
#### USAGE: awk -f generate-cprops-file.awk <path-to-fasta>
#### Output: cprops-formatted stdout
{
if ($0 ~ />/)
{
        if ( FNR != 1 )
		print old_id, new_id, clength;
	old_id = substr($1,2)
	new_id += 1
	clength = 0
}
else
         clength += length
}
END{
	print old_id, new_id, clength;
}
