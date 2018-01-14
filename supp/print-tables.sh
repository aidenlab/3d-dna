#!/bin/bash

pipeline=`cd "$( dirname $0)" && cd .. && pwd`

	bash ${pipeline}/supp/generate-table-1.sh
	
	bash ${pipeline}/supp/generate-table-s1.sh
	
	bash ${pipeline}/supp/generate-table-s2.sh
	
	bash ${pipeline}/supp/generate-table-s3.sh
	
	bash ${pipeline}/supp/generate-table-s4.sh
	
	bash ${pipeline}/supp/generate-table-s5.sh
