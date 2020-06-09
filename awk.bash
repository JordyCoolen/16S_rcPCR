FILE1=$1
awk 'BEGIN {FS="\t"; OFS="\t"} {print $5, $5, $2, $3, $4, $6}' $FILE1