#file=pair.txt
file=avedist2surface.txt
lines_in_file=`wc -l < $file`
lines_wanted=$(($lines_in_file/10))

shuf -n $lines_wanted $file
