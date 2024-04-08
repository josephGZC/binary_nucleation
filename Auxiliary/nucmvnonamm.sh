#!/bin/bash -v
# AUXILIARY SCRIPT
# DATE WRITTEN: SEPTEMBER 13, 2020
# DATE UPDATED: JUNE 6, 2021
# PURPOSE: rename NON-AMM minimization files used in NUCLEATION - NSRI PROJECT 

echo " "
echo " | AUXILIARY SCRIPT | RENAME EM NUCLEATION FILES | "
echo " "

# <=======================================================
# <==== I. DECLARE VARIABLES
# <=======================================================

molt1="NON"
molt2="AMM"

# <=======================================================
# <==== II. READ CURRENT FILES
# <=======================================================

ls *unitem_* > list1.txt
ls *unitem_* > list2.txt

# <=======================================================
# <==== III. SETTING DESIRED NAME
# <=======================================================

vim -c ":%s/unitem/unit_emm4/g" -c ":wq" list2.txt
vim -c ":%s/_box//g" -c ":wq" list2.txt
vim -c ":%s/ammonia/AMM/g" -c ":wq" list2.txt
vim -c ":%s/butanol/BUT/g" -c ":wq" list2.txt

# <=======================================================
# <==== IV. CHANGE NAME
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
echo " | TIME ELAPSED: $(($duration / 60)) MINUTE/S and $(($duration % 60)) SECOND/S |"
echo " "
