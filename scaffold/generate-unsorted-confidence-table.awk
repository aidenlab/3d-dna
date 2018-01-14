#### Description: Transform the density graph to a confidence graph.
#### Written by: OD
#### TODO: explore thresholds - raw coverage? misassembly upstream? repeats upstream? what else? BEGIN{SCORE_THRESH = 3}
{
counter[NR] = $0 # will need this info later

# initialize to zero values for top two partners
	left_best[$1]+=0
	left_secondbest[$1]+=0
	right_best[$1]+=0
	right_secondbest[$1]+=0
	left_best[$2]+=0
	left_secondbest[$2]+=0
	right_best[$2]+=0
	right_secondbest[$2]+=0

# keep track of the two top partners (on any side) for each of the contigs in the pair

if ($3 == 1)
{
	if ($4 >= left_best[$1])
	{
		left_secondbest[$1] = left_best[$1]
		left_best[$1] = $4
	}
	else if ($4 >= left_secondbest[$1])
		left_secondbest[$1] = $4

	if ($4 >= left_best[$2])
	{
		left_secondbest[$2] = left_best[$2]
		left_best[$2] = $4
	}
	else if ($4 >= left_secondbest[$2])
		left_secondbest[$2] = $4
}
else if ($3 == 2)
{
	if ($4 >= left_best[$1])
	{
		left_secondbest[$1] = left_best[$1]
		left_best[$1] = $4
	}
	else if ($4 >= left_secondbest[$1])
		left_secondbest[$1] = $4

	if ($4 >= right_best[$2])
	{
		right_secondbest[$2] = right_best[$2]
		right_best[$2] = $4
	}
	else if ($4 >= right_secondbest[$2])
		right_secondbest[$2] = $4
}
else if ($3 == 3)
{
	if ($4 >= right_best[$1])
	{
		right_secondbest[$1] = right_best[$1]
		right_best[$1] = $4
	}
	else if ($4 >= right_secondbest[$1])
		right_secondbest[$1] = $4

	if ($4 >= left_best[$2])
	{
		left_secondbest[$2] = left_best[$2]
		left_best[$2] = $4
	}
	else if ($4 >= left_secondbest[$2])
		left_secondbest[$2] = $4
}

else if ($3 == 4)
{
	if ($4 >= right_best[$1])
	{
		right_secondbest[$1] = right_best[$1]
		right_best[$1] = $4
	}
	else if ($4 >= right_secondbest[$1])
		right_secondbest[$1] = $4

	if ($4 >= right_best[$2])
	{
		right_secondbest[$2] = right_best[$2]
		right_best[$2] = $4
	}
	else if ($4 >= right_secondbest[$2])
		right_secondbest[$2] = $4
}

			
}
END{
for (n in counter) {
	split(counter[n],a," ")
	
	if (a[3] == 1)
	{
		data[1] = left_best[a[1]]
		data[2] = left_secondbest[a[1]]
		data[3] = left_best[a[2]]
		data[4] = left_secondbest[a[2]]
	}
	else if (a[3] == 2)
	{
		data[1] = left_best[a[1]]
		data[2] = left_secondbest[a[1]]
		data[3] = right_best[a[2]]
		data[4] = right_secondbest[a[2]]
	}
	else if (a[3] == 3)
	{
		data[1] = right_best[a[1]]
		data[2] = right_secondbest[a[1]]
		data[3] = left_best[a[2]]
		data[4] = left_secondbest[a[2]]
	}
	else if (a[3] == 4)
	{
		data[1] = right_best[a[1]]
		data[2] = right_secondbest[a[1]]
		data[3] = right_best[a[2]]
		data[4] = right_secondbest[a[2]]
	}

	asort(data) #sort the first and second best choices

	# correct for the double-counted [$1, $2] case
	if (data[3] == a[4])
		second = data[2]
	else
		second = data[3]

#TODO: 	introduce thresholding for inter vs intra? Can be done later as well

	if (second > 0)
	{
		conf = a[4]/second
	}
	else
#TODO: 	explore how to best deal with no-alternative links. Probably not too big of a problem for LIGer
		conf = 0
	print a[1], a[2], a[3], sprintf("%.20f",conf), sprintf("%.20f",a[4])
  	}
}


#
#
#
#
#
#
