/^[[:blank:]]*#/ {next} # ignore comments (lines starting with #)  
NF < 3 {next} # ignore lines which don’t have at least 3 columns  
$1 != prev {printf "\n" > "newf2000.txt"; prev=$1} # print blank line  
{print > "newf2000.txt"} # print the line
