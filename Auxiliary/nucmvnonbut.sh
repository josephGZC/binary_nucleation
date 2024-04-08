#!/bin/bash -v
# AUXILIARY SCRIPT
# DATE WRITTEN: SEPTEMBER 30, 2020
# DATE UPDATED: JULY 20, 2021
# PURPOSE: rename snapshot files used in NUCLEATION - NSRI PROJECT

echo " "
echo " | AUXILIARY SCRIPT | RENAME SNAPSHOT NUCLEATION FILES | "
echo " "

# <=======================================================
# <==== I. DECLARE VARIABLES
# <=======================================================

molt1="NON"
molt2="BUT"

# <=======================================================
# <==== II. READ CURRENT FILES
# <=======================================================

ls *.png > list1.txt
ls  ../../"Snapshots - Cross Section"/Cropped/*.jpg > list2.txt
sed 's/^........................................//' list2.txt > list3.txt

# <=======================================================
# <==== III. CHECK STATUS
# <=======================================================

check=$(head -1 list1.txt)
foo=$check 
for (( i=0; i<${#foo}; i++ )); do  #check if files already renamed
cchar=${foo:$i:1}
if [[ $cchar == "k" ]]; then
rm *list*
exit 0
fi
done

# <=======================================================
# <==== III. SETTING DESIRED NAME 
# <=======================================================

#vim -c ":%s/nb/NON_BUT/g" -c ":wq" list2.txt
#vim -c ":%s/ba/BUT_AMM/g" -c ":wq" list2.txt
#vim -c ":%s/na/NON_AMM/g" -c ":wq" list2.txt
#vim -c ":%s/wa/SOL_AMM/g" -c ":wq" list2.txt
#vim -c ":%s/wb/SOL_BUT/g" -c ":wq" list2.txt
#vim -c ":%s/wn/SOL_NON/g" -c ":wq" list2.txt
#vim -c ":%s/_VDW(0.1)/MER_cross/g" -c ":wq" list2.txt

#vim -c ":%s/MER_cross/MER/g" -c ":wq" list2.txt
#vim -c ":%s/_whole_/_/g" -c ":wq" list2.txt

#vim -c ":%s/_cross_/_/g" -c ":wq" list2.txt
#vim -c ":%s/MER/MER_cross/g" -c ":wq" list2.txt

#vim -c ":%s/NON/SOL/g" -c ":wq" list2.txt
#vim -c ":%s/BUT/AMM/g" -c ":wq" list2.txt
#vim -c ":%s/(/_/g" -c ":wq" list2.txt
#vim -c ":%s/)/MER/g" -c ":wq" list2.txt

#vim -c ":%s/_size/NS(/g" -c ":wq" list2.txt
#vim -c ":%s/time/${molt1}_${molt2}_/g" -c ":wq" list2.txt
#vim -c ":%s/_xaxis/)_xaxis/g" -c ":wq" list2.txt
#vim -c ":%s/_yaxis/)_yaxis/g" -c ":wq" list2.txt
#vim -c ":%s/_zaxis/)_zaxis/g" -c ":wq" list2.txt
#vim -c ":%s/.png/).png/g" -c ":wq" list2.txt
#vim -c ":%s/s).png/s.png/g" -c ":wq" list2.txt

file1="list1.txt"
file2="list3.txt"

# <=======================================================
# <==== IV. CHANGE NAME 
# <=======================================================

while read -r -u 3 file1 && read -r -u 4 file2; do
mv $file1 $file2
done 3<"$file1" 4<"$file2"

# ========================================================
# ------------------      END      -----------------------
# ========================================================

echo "...DONE"
echo " "

duration=$SECONDS
echo " | TIME ELAPSED: $(($duration / 60)) MINUTE/S and $(($duration % 60)) SECOND/S |"
echo " "
