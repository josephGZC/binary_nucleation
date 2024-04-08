#!/bin/bash -v
# AUXILIARY SCRIPT
# DATE WRITTEN: SEPTEMBER 13, 2020
# DATE UPDATED: JULY 20, 2021
# PURPOSE: rename SOL-BUT files used in NUCLEATION - NSRI PROJECT

echo " "
echo " | AUXILIARY SCRIPT | RENAME SIMULATION NUCLEATION FILES | "
echo " "

# A. UNIT - FIRST EM RUN

# <=======================================================
# <==== I. LIST PRESENT FILES
# <=======================================================

ls *unitem_* > list1.txt
ls *unitem_* > list2.txt

# <=======================================================
# <==== II. SETTING DESIRED NEW NAME
# <=======================================================

vim -c ":%s/unitem/unit_emm1/g" -c ":wq" list2.txt
vim -c ":%s/_box//g" -c ":wq" list2.txt
vim -c ":%s/butanol/BUT/g" -c ":wq" list2.txt
vim -c ":%s/water/SOL/g" -c ":wq" list2.txt

# <=======================================================
# <==== III. CHANGE NAME
# <=======================================================

file1="list1.txt"
file2="list2.txt"

while read -r -u 3 file1 && read -r -u 4 file2; do
mv $file1 $file2
done 3<"$file1" 4<"$file2"

# --------------      END PART 1      --------------------

# B. UNIT - SECOND TO FOURTH EM RUN

# <=======================================================
# <==== I. SET ITERATION
# <=======================================================

kip=1
while [ ${kip} -lt 4 ]; do

# <=======================================================
# <==== II. LIST PRESENT FILES
# <=======================================================

ls *em${kip}_* > list1.txt
ls *em${kip}_* > list2.txt

kap=$((${kip}+1))

# <=======================================================
# <==== III. SETTING DESIRED NEW NAME
# <=======================================================

vim -c ":%s/unitem${kip}/unit_emm${kap}/g" -c ":wq" list2.txt
vim -c ":%s/_box//g" -c ":wq" list2.txt
vim -c ":%s/butanol/BUT/g" -c ":wq" list2.txt
vim -c ":%s/water/SOL/g" -c ":wq" list2.txt

# <=======================================================
# <==== IV. CHANGE NAME
# <=======================================================

file1="list1.txt"
file2="list2.txt"

while read -r -u 3 file1 && read -r -u 4 file2; do
mv $file1 $file2
done 3<"$file1" 4<"$file2"

kip=$((${kip}+1))
done

# --------------      END PART 2      --------------------

# C. SYSTEM - FIRST EM RUN

# <=======================================================
# <==== I. LIST PRESENT FILES
# <=======================================================

ls *em.* > list1.txt    #for system_bsol_em.gro
ls *em.* > list2.txt

# <=======================================================
# <==== II. SETTING DESIRED NEW NAME
# <=======================================================

vim -c ":%s/system_bsol_em./system_emm1_SOL_BUT./g" -c ":wq" list2.txt

# <=======================================================
# <==== III. CHANGE NAME
# <=======================================================

file1="list1.txt"
file2="list2.txt"

while read -r -u 3 file1 && read -r -u 4 file2; do
mv $file1 $file2
done 3<"$file1" 4<"$file2"

# --------------      END PART 3      --------------------

# D. SYSTEM - SECOND TO FOURTH EM RUN

# <=======================================================
# <==== I. SET ITERATION
# <=======================================================

kip=1
while [ ${kip} -lt 4 ]; do

# <=======================================================
# <==== II. LIST PRESENT FILES                                                       
# <=======================================================

ls *_em${kip}* > list1.txt
ls *_em${kip}* > list2.txt

kap=$((${kip}+1))

# <=======================================================
# <==== III. SETTING DESIRED NEW NAME
# <=======================================================

vim -c ":%s/system_bsol_em${kip}/system_emm${kap}_SOL_BUT/g" -c ":wq" list2.txt

# <=======================================================                           
# <==== IV. CHANGE NAME
# <=======================================================

file1="list1.txt"
file2="list2.txt"

while read -r -u 3 file1 && read -r -u 4 file2; do
mv $file1 $file2
done 3<"$file1" 4<"$file2"

kip=$((${kip}+1))
done

# --------------      END PART 4      --------------------

# E. SYSTEM - NVT RUN

# <=======================================================
# <==== I. LIST PRESENT FILES
# <=======================================================

ls *_nvt.* > list1.txt    #for system_bsol_em.gro
ls *_nvt.* > list2.txt

# <=======================================================
# <==== II. SETTING DESIRED NEW NAME
# <=======================================================

vim -c ":%s/system_bsol_nvt/system_nvt1_SOL_BUT/g" -c ":wq" list2.txt

# <=======================================================
# <==== III. CHANGE NAME                                                             
# <=======================================================

file1="list1.txt"
file2="list2.txt"

while read -r -u 3 file1 && read -r -u 4 file2; do
mv $file1 $file2
done 3<"$file1" 4<"$file2"

# --------------      END PART 5      --------------------
                                                                                     
# F. SYSTEM - FIRST TO THIRD MD RUN

# <=======================================================
# <==== I. SET ITERATION
# <=======================================================

kip=1
while [ ${kip} -lt 4 ]; do

# <=======================================================
# <==== II. LIST PRESENT FILES
# <=======================================================

ls *_md${kip}.* > list1.txt
ls *_md${kip}.* > list2.txt

# <=======================================================
# <==== III. SETTING DESIRED NEW NAME
# <=======================================================

vim -c ":%s/system_bsol_md${kip}/system_mdd${kip}_SOL_BUT/g" -c ":wq" list2.txt

# <=======================================================
# <==== IV. CHANGE NAME
# <=======================================================

file1="list1.txt"
file2="list2.txt"

while read -r -u 3 file1 && read -r -u 4 file2; do
mv $file1 $file2
done 3<"$file1" 4<"$file2"

kip=$((${kip}+1))
done

# --------------      END PART 6      --------------------                           

# G. SYSTEM - FIRST TO THIRD MD RUN - CPT

# <=======================================================
# <==== I. SET ITERATION
# <=======================================================

kip=1
while [ ${kip} -lt 4 ]; do

# <=======================================================
# <==== II. LIST PRESENT FILES
# <=======================================================

ls *_md${kip}_* > list1.txt
ls *_md${kip}_* > list2.txt

# <=======================================================
# <==== III. SETTING DESIRED NEW NAME
# <=======================================================

vim -c ":%s/system_bsol_md${kip}/system_mdd${kip}_SOL_BUT/g" -c ":wq" list2.txt

# <=======================================================
# <==== IV. CHANGE NAME
# <=======================================================

file1="list1.txt"
file2="list2.txt"

while read -r -u 3 file1 && read -r -u 4 file2; do
mv $file1 $file2
done 3<"$file1" 4<"$file2"

kip=$((${kip}+1))
done

# --------------      END PART 7      --------------------

# H. SYSTEM - FIRST TO THIRD MD RUN - TRR
                                                                                     
# <=======================================================
# <==== I. SET ITERATION
# <=======================================================

kip=2
while [ ${kip} -lt 4 ]; do

# <=======================================================
# <==== II. LIST PRESENT FILES
# <=======================================================

ls *_md1-* > list1.txt
ls *_md1-* > list2.txt

vim -c ":%s/system_bsol_md1-${kip}/system_mdd1-${kip}_SOL_BUT/g" -c ":wq" list2.txt

file1="list1.txt"
file2="list2.txt"

while read -r -u 3 file1 && read -r -u 4 file2; do
mv $file1 $file2
done 3<"$file1" 4<"$file2"

kip=$((${kip}+1))
done

# --------------      END PART 8      --------------------

# ========================================================
# ------------------      END      -----------------------
# ========================================================  

echo "...DONE"
echo " "

duration=$SECONDS
echo " | TIME ELAPSED: $(($duration / 60)) MINUTE/S and $(($duration % 60)) SECOND/S |"
echo " "
