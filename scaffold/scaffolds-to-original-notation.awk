#### Description: Helper script to replace scaffold pool.
#### Written by: OD
function abs(value)
{
  return (value<0?-value:value);
}
function reverse(seq)
{
  split (seq,s)
  tmpstring = -s[length(s)]
  for (i = length(s)-1; i >= 1; i--)
	tmpstring = sprintf("%s %s",tmpstring, -s[i])
  return tmpstring;
}
{
	if (FILENAME==ARGV[1])
	{
	for (n = 1; n <= NF; n++)
		{
		scaffold[abs($n)] = NR
		position[abs($n)] = ($n/abs($n))*n
		}
	max_nf[NR] = NF
	max_nr = NR
	next
	}
}
{
	if (FNR in scaffold)
	{	
		if (position[FNR] > 0) 
			newscaffold[scaffold[FNR], abs(position[FNR])] = $0
		else
			newscaffold[scaffold[FNR], abs(position[FNR])] = reverse($0)

	}
	else
	{
		max_nr += 1
		newscaffold[max_nr, 1] = $0
		max_nf[max_nr] = 1
	}

}
END{
	for (n = 1; n <= max_nr; n++)
		{	
			tmpstring = newscaffold[n, 1]
			if (max_nf[n] > 1)
				{
				for (k = 2; k <= max_nf[n]; k++)
					tmpstring = sprintf("%s %s", tmpstring, newscaffold[n, k])
				}
			print tmpstring
		}	
   }
