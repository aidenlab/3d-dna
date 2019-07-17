#### Description: path cover from confidence file.
#### USAGE: awk -f confidence-to-assembly.awk <path-to-confidence-file>
#### Output: asm-formatted stdout
#### Written by: OD
function abs(value)
{
  return (value<0?-value:value);
}
function reverse(seq)
{
  split (seq,s)
  tmpstring = -s[length(s)]
  for (i = length(s)-1; i>=1; i--)
	tmpstring = sprintf("%s %s",tmpstring, -s[i])
  return tmpstring;
}
BEGIN{CONFTHRES=1}
{
if (NR == 1)
{	if($4 > CONFTHRES)
	{
		if ($3 == 1)
			sequence[1] = sprintf("%s %s", -$2, $1)
		else if	($3 == 2)
			sequence[1] = sprintf("%s %s", $2, $1)
		else if	($3 == 3)
			sequence[1] = sprintf("%s %s", $1, $2)
		else if ($3 == 4)
			sequence[1] = sprintf("%s %s", $1, -$2)
		scaf[$1] = 1
		scaf[$2] = 1
	}
	else {exit 0}
}
else
{
	if ($4 <= CONFTHRES) {exit 0}
	else 	{
		pos1 = scaf[$1]
		pos2 = scaf[$2]
		}
		if ((pos1 == "") && (pos2 == ""))
			{
			newindex = length(sequence) + 1;
			while (newindex in sequence) {newindex++}
			if ($3 == 1)
				sequence[newindex] = sprintf("%s %s", -$2, $1)
			else if	($3 == 2)
				sequence[newindex] = sprintf("%s %s", $2, $1)
			else if	($3 == 3)
				sequence[newindex] = sprintf("%s %s", $1, $2)
			else if ($3 == 4)
				sequence[newindex] = sprintf("%s %s", $1, -$2)
			scaf[$1] = newindex
			scaf[$2] = newindex	
			}
		else if ((pos1 !="") && (pos2 !=""))
		{
			if (pos1 != pos2)
			{
				split(sequence[pos1],a)
				split(sequence[pos2],b)
				first1 = a[1]
				last1 = a[length(a)]
				first2 = b[1]
				last2 = b[length(b)]
				if ($3 == 1)
				{
					skip = 0
					if ((first1 == $1) && (last2 == -$2))
						sequence[pos1] = sprintf("%s %s", sequence[pos2], sequence[pos1])			
					else if ((last1 == -$1) && (first2 == $2))
						sequence[pos1] = sprintf("%s %s", sequence[pos1], sequence[pos2])			
					else if ((last1 == -$1) && (last2 == -$2))
						sequence[pos1] = sprintf("%s %s", sequence[pos1], reverse(sequence[pos2]))		
					else if ((first1 == $1) && (first2 == $2))
						sequence[pos1] = sprintf("%s %s", reverse(sequence[pos2]), sequence[pos1])		
					else skip = 1
				}
					
				else if ($3 == 2)
				{
					skip = 0
					if ((first1 == $1) && (last2 == $2))
						sequence[pos1] = sprintf("%s %s", sequence[pos2], sequence[pos1])
					else if ((last1 == -$1) && (first2 == -$2))
						sequence[pos1] = sprintf("%s %s", sequence[pos1], sequence[pos2])
					else if ((last1 == -$1) && (last2 == $2))
						sequence[pos1] = sprintf("%s %s", sequence[pos1], reverse(sequence[pos2]))
					else if ((first1 == $1) && (first2 == -$2))
						sequence[pos1] = sprintf("%s %s", reverse(sequence[pos2]), sequence[pos1])
					else skip = 1
				}
				
				else if ($3 == 3)
				{
					skip = 0
					if ((last1 == $1) && (first2 == $2))
						sequence[pos1] = sprintf("%s %s", sequence[pos1], sequence[pos2])
					else if ((first1 == -$1) && (last2 == -$2))
						sequence[pos1] = sprintf("%s %s", sequence[pos2], sequence[pos1])
					else if ((last1 == $1) && (last2 == -$2))
						sequence[pos1] = sprintf("%s %s", sequence[pos1], reverse(sequence[pos2]))
					else if ((first1 == -$1) && (first2 == $2))
						sequence[pos1] = sprintf("%s %s", reverse(sequence[pos2]), sequence[pos1])
					else skip = 1
				}
				
				else if ($3 == 4)
				{
					skip = 0
					if ((last1 == $1) && (first2 == -$2))
						sequence[pos1] = sprintf("%s %s", sequence[pos1], sequence[pos2])
					else if ((first1 == -$1) && (last2 == $2))
						sequence[pos1] = sprintf("%s %s", sequence[pos2], sequence[pos1])
					else if ((last1 == $1) && (last2 == $2))
						sequence[pos1] = sprintf("%s %s", sequence[pos1], reverse(sequence[pos2]))
					else if ((first1 == -$1) && (first2 == -$2))
						sequence[pos1] = sprintf("%s %s", reverse(sequence[pos2]), sequence[pos1])
					else skip = 1
				}

				if (skip != 1)
				{	
					for (i in b)
						scaf[abs(b[i])] = pos1
					delete sequence[pos2]
				}

			}
		}
		else if ((pos1 != "") && (pos2 == ""))
		{			
			split(sequence[pos1],a)
			first1 = a[1]
			last1 = a[length(a)]
			skip = 0
			if ($3 == 1)
				if (first1 == $1)
					sequence[pos1] = sprintf("%s %s", -$2, sequence[pos1])
				else if (last1 == -$1)
					sequence[pos1] = sprintf("%s %s", sequence[pos1], $2)
				else skip = 1
			else if ($3 == 2)
				if (first1 == $1)
					sequence[pos1] = sprintf("%s %s", $2, sequence[pos1])
				else if (last1 == -$1)
					sequence[pos1] = sprintf("%s %s", sequence[pos1], -$2)
				else skip = 1
			else if ($3 == 3)
				if (last1 == $1)
					sequence[pos1] = sprintf("%s %s", sequence[pos1], $2)
				else if (first1 == -$1)
					sequence[pos1] = sprintf("%s %s", -$2, sequence[pos1])
				else skip = 1
			else if ($3 == 4)
				if (last1 == $1)
					sequence[pos1] = sprintf("%s %s", sequence[pos1], -$2)
				else if (first1 == -$1)
					sequence[pos1] = sprintf("%s %s", $2, sequence[pos1])
				else skip = 1
			if (skip != 1)
				scaf[$2] = pos1
		}
		else if ((pos1 == "") && (pos2 != ""))
		{
			split(sequence[pos2],b)
			first2 = b[1]
			last2 = b[length(b)]
			skip = 0
			if ($3 == 1)
				if (last2 == -$2)
					sequence[pos2] = sprintf("%s %s", sequence[pos2], $1)
				else if (first2 == $2)
					sequence[pos2] = sprintf("%s %s", -$1, sequence[pos2])
				else skip = 1
			else if ($3 == 2)
				if (last2 == $2)
					sequence[pos2] = sprintf("%s %s", sequence[pos2], $1)
				else if (first2 == -$2)
					sequence[pos2] = sprintf("%s %s", -$1, sequence[pos2])
				else skip = 1
			else if ($3 == 3)
				if (first2 == $2)
					sequence[pos2] = sprintf("%s %s", $1, sequence[pos2])
				else if (last2 == -$2)
					sequence[pos2] = sprintf("%s %s", sequence[pos2], -$1)
				else skip = 1
			else if ($3 == 4)
				if (first2 == -$2)
					sequence[pos2] = sprintf("%s %s", $1, sequence[pos2])
				else if (last2 == $2)
					sequence[pos2] = sprintf("%s %s", sequence[pos2], -$1)
				else skip = 1
			if (skip != 1)
				scaf[$1] = pos2
		}
}
}
END{for (n in sequence) 
	{
		print sequence[n]
	}
   }
