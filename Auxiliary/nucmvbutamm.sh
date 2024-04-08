#!/bin/bash -v
# AUXILIARY SCRIPT
# DATE WRITTEN: SEPTEMBER 13, 2020 
# DATE UPDATED: JULY 20, 2021
# PURPOSE: rename ammbut system name to BUT_AMM in NUCLEATION - NSRI PROJECT

echo " "
echo " | AUXILIARY SCRIPT | RENAME NUCLEATION SIMULATION FILES | "
echo " "

# <=======================================================
# <==== I. LIST PRESENT FILES
# <=======================================================

ls *system_* > list1.txt
ls *system_* > list2.txt

# <=======================================================
# <==== II. SETTING DESIRED NEW NAME
# <=======================================================

vim -c ":%s/ammbut/BUT_AMM/g" -c ":wq" list2.txt

# <=======================================================
# <==== III. CHANGE NAME
# <=======================================================

file1="list1.txt"
file2="list2.txt"

while read -r -u 3 file1 && read -r -u 4 file2; do
mv $file1 $file2
done 3<"$file1" 4<"$file2"

# ========================================================
# ------------------      END      -----------------------
# ========================================================

echo "...DONE"
echo " "

duration=$SECONDS
echo " | TIME ELAPSED: $(($duration / 60)) MINUTE/S and $(($duration % 60)) SECOND/S | "
echo " " 
