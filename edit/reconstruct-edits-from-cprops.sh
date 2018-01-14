#!/bin/bash

pipeline=`cd "$( dirname $0)" && cd .. && pwd`

orig_cprops=$1
current_cprops=$2

awk 'BEGIN{OFS="\t"; print "chr1", "x1", "x2", "chr2", "y1", "y2", "color", "id", "X1", "X2", "Y1", "Y2"}$1~/:::debris/{print $1, 0, $3, $1, 0, $3, "0,0,0", "debris", 0, $3, 0, $3}' ${current_cprops} | awk -f ${pipeline}/lift/lift-input-annotations-to-asm-annotations.awk ${current_cprops} <(awk '{print $2}' ${current_cprops}) - | awk -f ${pipeline}/lift/lift-asm-annotations-to-input-annotations.awk ${orig_cprops} <(awk '{print $2}' ${orig_cprops}) - 
