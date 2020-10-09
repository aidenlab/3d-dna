## Helper script to parse a fasta for components. Should maybe get rid of this in favor of construct-fasta-from-asm.sh
FILENAME==ARGV[1]{include[$1];next}
$0~/>/{test=0; if(substr($1,2) in include){test=1}; if(except){test=!test}}
test{print}
