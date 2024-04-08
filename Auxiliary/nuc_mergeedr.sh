#!/bin/bash
# AUXILIARY SCRIPT
# DATE WRITTEN: AUGUST 21, 2020
# DATE UPDATED: AUGUST 30, 2021
# PURPOSE: output energy report

echo " "
echo " | AUXILIARY SCRIPT | OUTPUT ENERGY REPORT | "
echo " "

# ========================================================
# ------------------    REMINDERS    ---------------------
# ========================================================

# 1. only works for production files labelled as: system_mdd*
#
# 2. run script in csrc as
#    nohup ./nuc_mergeedr.sh > nohup.out &
#
# 3. the following files are generated by this script:
#    a. ${molec1}_${molec2}_SYSREPORT.edr          (gmx enecov) 
#    b. ${molec1}_${molec2}_SYSREPORT1.log         (gmx energy) [PE, KE, TOTAL]
#    c. ${molec1}_${molec2}_SYSREPORT2.log         (gmx energy) [TEMP]
#    d. ${molec1}_${molec2}_ALL-ENERGY1.xvg & .png (gmx energy) [PE, KE, TOTAL]
#    e. ${molec1}_${molec2}_ALL-ENERGY2.xvg & .png (gmx energy) [TEMP]

# ========================================================
# -------------------    SECTIONS    ---------------------
# ========================================================

# I.   DETERMINE SYSTEM
#
# II.  CORRECT SYNTAX
#
# III. REGISTER FILE ID
#
# IV.   REGISTER ENERGY FILENAMES
#
# V.    TALLY ALL ENERGY FILES
#
# VI.   CONCATENATE ALL EDR FILES
#
# VII.  OUTPUT DESIRED ENERGY REPORTS
#
# VIII. SETUP GRAPH SETTINGS
#
# IX.   CONVERT XVG TO PNG
#
# X.    TRANSFER REPORTS IN DESIGNATED DIRECTORY

# <=======================================================
# <==== I. DETERMINE SYSTEM
# <=======================================================

ls system_mdd*.edr > list1.txt

detnme=$( tail -1 list1.txt )
molec1=${detnme:12:3}
molec2=${detnme:16:3}

# <=======================================================
# <==== II. CORRECT SYNTAX
# <=======================================================

[[ ! -f "system_mdd1_${molec1}_${molec2}.part0001.edr" ]] &&
mv system_mdd1_${molec1}_${molec2}.edr system_mdd1_${molec1}_${molec2}.part0001.edr 

# <=======================================================
# <==== III. REGISTER FILE ID
# <=======================================================

compvl=1
while read -r row
do
compvt=${row:10:1}
if [[ $compvt -gt $compvl ]]; then
compvl=$compvt            # highest number 
fi
done < list1.txt
rm list1.txt

# <=======================================================
# <==== IV. REGISTER ENERGY FILENAMES
# <=======================================================

currvl=0
while [[ $compvl -ne $currvl ]] ; do

currvl=$((currvl + 1))
checkc="000${currvl}"
alignc=$(echo $checkc | rev | cut -b 1-4 | rev)

prodnc="system_mdd${currvl}_${molec1}_${molec2}.part${alignc}.edr"

# <=======================================================
# <==== V. TALLY ALL ENERGY FILES
# <=======================================================

if [[ prodnc == 1 ]]; then
basnam="$prodnc"
else
basnam="$basnam $prodnc"
fi
done

# <=======================================================
# <==== VI. CONCATENATE ALL EDR FILES
# <=======================================================

gmx eneconv -f ${basnam} -o ${molec1}_${molec2}_SYSREPORT.edr

[[ ! -f "${molec1}_${molec2}_SYSREPORT.edr" ]] &&
echo "ERROR [1]: ${molec1}_${molec2}_SYSREPORT.edr DOES NOT EXIST"  && exit

# <=======================================================
# <==== VII. OUTPUT DESIRED ENERGY REPORTS
# <=======================================================

gmx energy -f ${molec1}_${molec2}_SYSREPORT.edr -o ${molec1}_${molec2}_ALL-ENERGY1.xvg << EOF >> ${molec1}_${molec2}_SYSREPORT1.log
Pot
Kin
Tot
EOF

gmx energy -f ${molec1}_${molec2}_SYSREPORT.edr -o ${molec1}_${molec2}_ALL-ENERGY2.xvg << EOF >> ${molec1}_${molec2}_SYSREPORT2.log
Tem
EOF

# <=======================================================
# <==== VIII. SETUP GRAPH SETTINGS
# <=======================================================

for file in *ALL-ENERGY*.xvg ; do
docu=$( echo $file | rev | cut -c 5- | rev )
mark=$( echo $file | rev | cut -c 5  | rev )

SED_COM="/title/ a "
SED_COM+="@ title font 6 \n"
SED_COM+="@ subtitle \"${molec1}-${molec2}\" \n"
SED_COM+="@ subtitle font 6 \n" 
SED_COM+="@ xaxis label font 6 ; xaxis ticklabel font 6 \n"
SED_COM+="@ yaxis label font 6 ; yaxis ticklabel font 6 \n"
SED_COM+="@ xaxis tick minor linewidth 2 ; xaxis tick major linewidth 2 \n"         
SED_COM+="@ yaxis tick minor linewidth 2 ; yaxis tick major linewidth 2 \n"
SED_COM+="@ frame linewidth 2"
sed -i "${SED_COM}" $file

SED_COM="/length/ a "
SED_COM+="@ legend 0.187, 0.81 \n"
SED_COM+="@ legend font 6 \n"
SED_COM+="@ legend vgap 2 \n"
SED_COM+="@ AUTOTICKS \n"
SED_COM+="@ autoscale onread none \n"
SED_COM+="@ hardcopy device \"PNG\" \n"
SED_COM+="@ print to \"${docu}.png\""
sed -i "${SED_COM}" $file

[[ $mark == 1 ]] &&
sed -i "/AUTOTICKS/ i @ world 0,-300000,6000,500000" $file ||
sed -i "/AUTOTICKS/ i @ world 0,200,6000,700" $file
sed -i '/view /d' $file
sed -i '/legend 0.78/d' $file
done

# <=======================================================
# <==== VI. CONVERT XVG TO PNG 
# <=======================================================

grace -nxy ${molec1}_${molec2}_ALL-ENERGY1.xvg -hardcopy 
grace -nxy ${molec1}_${molec2}_ALL-ENERGY2.xvg -hardcopy 

# <=======================================================
# <==== VII. TRANFER REPORTS IN DESIGNATED DIRECTORY
# <=======================================================

[[ ! -d ./energy/xvg ]] &&
mkdir ./energy
mv *_SYSREPORT* ./energy

[[ ! -d ./energy/xvg ]] &&
mkdir ./energy/xvg
mv *_ALL-ENERGY*.xvg ./energy/xvg

[[ ! -d ./energy/png ]] &&
mkdir ./energy/png
mv *_ALL-ENERGY*.png ./energy/png

echo "...ENERGY DIRECTORY CONTENTS"

echo " "
echo "...xvg"
ls -ltr ./energy/xvg | awk '{print $6, $7, $8, $9}' 

echo " "
echo "...png"
ls -ltr ./energy/png | awk '{print $6, $7, $8, $9}' 

# ========================================================
# ------------------      END      -----------------------
# ========================================================

echo " "
echo "...DONE"
echo " "

duration=$SECONDS
echo " | TIME ELAPSED: $(($duration / 60)) MINUTE/S and $(($duration % 60)) SECOND/S |"
echo " "