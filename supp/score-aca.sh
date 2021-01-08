#!/bin/bash

#### Description: Script to calculate scores from the aggregate chromosome analysis (ACA) maps.
#### Usage: bash ./score-aca.sh <path_to_aca_hic_file>.
#### Input: ACA Hi-C contact map in hic format.
#### Output: stdout.
#### Parameters: none [TODO: resolution at which the calculation is performed can be made into a parameter].
#### Dependencies: 3D-DNA scripts.
#### Written by: Olga Dudchenko, version date 08/04/2020

## 3D-DNA
pipeline=`cd "$( dirname $0)" && cd .. && pwd`

hic=$1

res=50000
#res=100000

o=$(bash ${pipeline}/visualize/juicebox_tools.sh dump observed KR ${hic} assembly:500000:500000 assembly:1500000:1500000 BP ${res} | awk 'END{printf "%f", $NF}')

e1=$(bash ${pipeline}/visualize/juicebox_tools.sh dump observed KR ${hic} assembly:$((500000+res)):$((500000+res)) assembly:$((1500000+res)):$((1500000+res)) BP ${res} | awk 'END{printf "%f", $NF}')

e2=$(bash ${pipeline}/visualize/juicebox_tools.sh dump observed KR ${hic} assembly:500000:500000 assembly:$((1500000+res)):$((1500000+res)) BP ${res} | awk 'END{printf "%f", $NF}')

Cscore=$(printf '%.*f\n' 3 $(bc -l <<< "2*$o/($e1+$e2)"))

o=$(bash ${pipeline}/visualize/juicebox_tools.sh dump observed KR ${hic} assembly:0:0 assembly:$((1000000-res)):$((1000000)) BP ${res} | awk '{c+=$NF}END{printf "%f", c}')

#o=$(bash ${pipeline}/visualize/juicebox_tools.sh dump observed KR ${hic} assembly:0:0 assembly:$((1000000)):$((1000000)) BP ${res} | awk 'END{printf "%f", $NF}')

e=$(bash ${pipeline}/visualize/juicebox_tools.sh dump observed KR ${hic} assembly:${res}:${res} assembly:$((1000000-2*res)):$((1000000+res)) BP ${res} | awk '{c+=$NF}END{printf "%f", c}')

#e=$(bash /${pipeline}/visualize/juicebox_tools.sh dump observed KR ${hic} assembly:${res}:${res} assembly:$((1000000)):$((1000000+res)) BP ${res} | awk '{c+=$NF}END{printf "%f", c}')

#e2=$(bash ${pipeline}/visualize/juicebox_tools.sh dump observed KR ${hic} assembly:${res}:${res} assembly:$((2000000-2*res)):$((2000000-2*res)) BP ${res} | awk 'END{printf "%f", $NF}')

Tscore=$(printf '%.*f\n' 3 $(bc -l <<< "2*$o/$e"))

o=$(bash ${pipeline}/visualize/juicebox_tools.sh dump observed KR ${hic} assembly:$((1000000-res)):$((1000000-res)) assembly:$((6*res)):$((1000000-7*res)) BP ${res} | awk '{c+=$NF}END{printf "%f", c}')

e=$(bash ${pipeline}/visualize/juicebox_tools.sh dump observed KR ${hic} assembly:1000000:1000000 assembly:$((6*res)):$((1000000-7*res)) BP ${res} | awk '{c+=$NF}END{printf "%f", c}')

CTscore=$(printf '%.*f\n' 3 $(bc -l <<< "$o/$e"))

o=$(bash ${pipeline}/visualize/juicebox_tools.sh dump observed KR ${hic} assembly:0:$((1000000-res)) assembly:0:$((1000000-res)) BP ${res} | awk -v res=${res} '($1+$2==1000000-res)&&$1>0*res&&$1<(500000-1*res){c+=$NF}END{printf "%f", c}')

e=$(bash ${pipeline}/visualize/juicebox_tools.sh dump observed KR ${hic} assembly:0:$((1000000-res)) assembly:0:$((1000000-res)) BP ${res} | awk -v res=${res} '($1+$2==1000000-3*res)&&$1>-1*res&&$1<(500000-2*res){c+=$NF}END{printf "%f", c}')

Fscore=$(printf '%.*f\n' 3 $(bc -l <<< "$o/$e"))

printf '%s\t%s\t%s\t%s\n' ${CTscore} ${Fscore} ${Cscore} ${Tscore}
