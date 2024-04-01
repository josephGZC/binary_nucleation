#!/bin/bash
# PART4.A OF ANALYSIS SHELL SCRIPTS
# DATE WRITTEN: AUGUST 27, 2020
# DATE UPDATED: FEBRUARY 9, 2022
# PURPOSE: generate .gro cross section structure files
SECONDS=0

echo " "
echo "         BRANCH TASK: cut.sh"
echo " "

# <=======================================================
# <==== I. LOOP ALL STRUCTURE FILES
# <=======================================================

for i in {1..5..1}
do

for j in {1..3..1}
do

rmd=5000
mxv=$(expr $i \* 10)
rmc=$(expr $rmd + $mxv + $j)

rawfil="fort.${rmc}"

# <=======================================================
# <==== II. DECLARE VARIABLES
# <=======================================================

if [[ -f $rawfil ]]; then
molety=$(awk '{if(NR==3) print $9}' ./graphs/gnu/PLOTfrequency.gnu)
timest=$(awk '{if(NR==1) print $4}' $rawfil | cut -c 1-5) 
mersze=$(awk '{if(NR==1) print $7}' $rawfil)
fileid="${molety}_${timest}NS_${mersze}MER_SNAPT"
groext="${fileid}.gro"
clean1="${fileid}_clean1.gro"
clean2="${fileid}_clean2.gro"
xaxisf="${fileid}_xaxis.gro"
yaxisf="${fileid}_yaxis.gro"
zaxisf="${fileid}_zaxis.gro"

# <=======================================================
# <==== III. RENAME STRUCTURE FILE 
# <=======================================================

cp $rawfil $groext 	
echo "   $rawfil to $groext"
echo " "

# <=======================================================
# <==== IV. CLEAN GRO FILE FOR CUTTING 
# <=======================================================

SED_COM="s/^.\{5\}//;"                        #removes first 5 characters                                                    
SED_COM+="s/\(.\{3\}\).\{12\}/\1   /;"        #replaces 3rd to 12 character to white space
SED_COM+="p"                                  #print

sed '/[=]/d;/[A-Z]/!d' $groext | sed -n "$SED_COM" > $clean1  #delete lines with equal sign, and lines 
                                                              #without capital letter
# <=======================================================
# <==== V. GENERATE SLICED CLUSTERS 
# <=======================================================

( sed -n '2{p;q;}' $rawfil && cat $clean1 ) > temporary && mv temporary $clean1     # inserting atom count
gfortran -o3 crosscutgro.f -o crosscutgro
./crosscutgro < $clean1      						   # generating cross section
cp "fort.60" $xaxisf	     						   
cp "fort.61" $yaxisf	     						   
cp "fort.62" $zaxisf	     						   
rm "fort.60" "fort.61" "fort.62"

fi

done
done

# =========================================================
# ------------------------- END ---------------------------
# =========================================================

duration=$SECONDS
echo "         TIME ELAPSED: $(($duration / 60)) MINUTE/S $(($duration % 60)) SECOND/S"
echo " "
