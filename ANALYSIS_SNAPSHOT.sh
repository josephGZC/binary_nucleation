#!/bin/bash
# PART4 OF ANALYSIS SCRIPTS 
# DATE WRITTEN: AUGUST 27, 2020                                             
# DATE UPDATED: JULY 20, 2021
# PURPOSE: prepare vmd rendering script
SECONDS=0

echo " "
echo " | COMMENCE: PART4 OF ANALYSIS SCRIPTS | PREPARE VMD RENDERING SCRIPTS | "
echo " "

# ========================================================
# ------------------    REMINDERS    ---------------------
# ========================================================
                                     
# 1. the sequence of analysis scripts are as follow:
#    [1] ./ANALYSIS_CLEAN.sh
#    [2] ./ANALYSIS_EXAMINE.sh   
#    [3] ./ANALYSIS_PLOT.sh
#    [4] ./ANALYSIS_SNAPSHOT.sh    (*)
#
# 2. further breakdown of the scripts are as follow:
#    [1] ./ANALYSIS_CLEAN.sh 
#
#    [2] ./ANALYSIS_EXAMINE.sh
#         (A) cluster_nm.f
#
#    [3] ./ANALYSIS_PLOT.sh 
#         (A) graphing.f
#         (B) PLOTRDP1A.gnu
#         (C) PLOTRDP1B.gnu
#         (D) PLOTRDP2A.gnu
#         (E) PLOTRDP2B.gnu
#         (F) PLOTRDP3A.gnu
#         (G) PLOTRDP3B.gnu
#         (H) PLOTRDF1A.gnu
#         (I) PLOTRDF1B.gnu
#         (J) PLOTfrequency.gnu
#         (K) PLOTmolfraction.gnu
#         (L) PLOTmass.gnu
#         (M) PLOTmassfracmax.gnu
#         (N) PLOTmassfracmin.gnu
#         (O) PLOTradius.gnu
#
#    [4] ./ANALYSIS_SNAPSHOT.sh (adjust scale variable sclset)
#         (A) cut.sh  
#                (a) crosscutgro.f
#         (B) generates renindivwnba.tcl
#             - activate using VMD tkconsole
#         (C) generates renbatchwnba.tcl
#             - activate using VMD tkconsole
#
#        Manual Command:
#        $ cd snapshots
#        $ vmd -dispdev none -e renbatchwnba.tcl
#        $ mogrify -antialias -density 300 -format png -quality 100 -colorspace RGB *.tg
#
# 3. the following files are required to be in the same directory
#    a. fort.10
#    b. fort.5011-13
#       fort.5021-23
#       fort.5031-33
#       fort.5041-43
#       fort.5051-53
#
# 4. declare important variables in section I, for:
#    a. scale of image rendered
# 
# 5. if script is run in:
#    a. CSRC
#       nohup "./ANALYSIS_EXAMINE.sh" > out.log &
#    b. ASTI
#       sbatch ANALYSIS_EXAMINE.sub

# ========================================================
# -------------------    SECTIONS    ---------------------
# ========================================================

# I.    DECLARE VARIABLES
#
# II.   GENERATE SNAPSHOT-COMPATIBLE FILES
#
# III.  CREATE WORKING DIRECTORY
#
# IV.   DETERMINE PAIR 
#
# V.    LIST STRUCTURE FILES
#
# VI.   GENERATE VMD SCRIPT HEADING
#
# VII.  WRITE STRUCTURE FILES IN VMD SCRIPT
#
# VIII. SETTING UP SCRIPT LOOP
#
# IX.   SETTING UP SCRIPT COMMANDS
#
# X.    GENERATE SCRIPT FOR INDIVIDUAL RENDERING
#
# XI.   REMOVE EXCESS FILES

# <=======================================================
# <==== I. DECLARE VARIABLES
# <=======================================================

echo "...(1/11) DECLARE VARIABLES"
echo " "

sclset=0.022

# <=======================================================
# <==== II. GENERATE SNAPSHOT-COMPATIBLE FILES
# <=======================================================

echo "...(2/11) GENERATE SNAPSHOT-COMPATIBLE FILES"
./cut.sh

# <=======================================================
# <==== III. CREATE WORKING DIRECTORY
# <=======================================================

echo "...(3/11) CREATE WORKING DIRECTORY"
echo " "
if [[ ! -d "snapshots/" ]]; then
mkdir snapshots
fi
mv *SNAPT* snapshots/
rm snapshots/*clean*

if [[ -f snapshots/*.eps ]]; then
rm snapshots/*.eps
fi

if [[ -f snapshots/renbatchwnba.tcl ]]; then
rm snapshots/renbatchwnba.tcl
fi

# <=======================================================
# <==== IV. DETERMINE PAIR 
# <=======================================================

echo "...(4/11) DETERMINE PAIR"
echo " "
typem1=$(awk '{if(NR==3) print $9}' ./graphs/gnu/PLOTfrequency.gnu | cut -c 1-3)
if [[ $typem1 == "SOL" ]]; then
molec1="SOL"	
color1="blue2"
elif [[ $typem1 == "NON" ]]; then
molec1="NON"	
color1="green3"
elif [[ $typem1 == "BUT" ]]; then
molec1="BUTA"	
color1="yellow"
elif [[ $typem1 == "AMM" ]]; then
molec1="NH3"	
color1="red2"
elif [[ $typem1 == "MET" ]]; then
molec1="MET"	
color1="orange2"
elif [[ $typem1 == "ACT" ]]; then
molec1="ACT"	
color1="pink"
elif [[ $typem1 == "OCT" ]]; then
molec1="OCTA"	
color1="pink"
fi

typem2=$(awk '{if(NR==3) print $9}' ./graphs/gnu/PLOTfrequency.gnu | cut -c 5-7)
if [[ $typem2 == "SOL" ]]; then
molec2="SOL"	
color2="blue2"
elif [[ $typem2 == "NON" ]]; then
molec2="NON"	
color2="green3"
elif [[ $typem2 == "BUT" ]]; then
molec2="BUTA"	
color2="yellow"
elif [[ $typem2 == "AMM" ]]; then
molec2="NH3"	
color2="red2"
elif [[ $typem2 == "MET" ]]; then
molec2="MET"	
color2="orange2"
elif [[ $typem2 == "ACT" ]]; then
molec2="ACT"	
color2="pink"
elif [[ $typem2 == "OCT" ]]; then
molec2="OCTA"	
color2="pink"
fi

# <=======================================================
# <==== V. LIST STRUCTURE FILES
# <=======================================================

echo "...(5/11) LIST STRUCTURE FILES"
echo " "
ls snapshots/*SNAPT* > snapshots/temp1.txt
sed 's/....$//; s/^..........//' snapshots/temp1.txt > snapshots/temp2.txt
rm snapshots/temp1.txt

# <=======================================================
# <==== VI. GENERATE VMD SCRIPT HEADING
# <=======================================================

echo "...(6/11) GENERATE VMD SCRIPT HEADING"
echo " "
cat <<EOT >> snapshots/renbatchwnba.tcl
# WNBA CLUSTER SNAPSHOT BATCH RENDERING SCRIPT 08/28/20
# CLUSTER SNAPSHOTS ARE RENDERED INTO TGA FILES WITH RES 2000 x 2000

# I. STRUCTURE FILES
EOT

# <=======================================================
# <==== VII. WRITE STRUCTURE FILES IN VMD SCRIPT
# <=======================================================

echo "...(7/11) WRITE STRUCTURE FILES IN VMD SCRIPT"
echo " "
cn=0
while read -r row 	
do
cn=$((cn+1))
echo "set fileid"${cn}" "${row}"" >> snapshots/renbatchwnba.tcl
done < snapshots/temp2.txt

# <=======================================================
# <==== VIII. SETTING UP SCRIPT LOOP
# <=======================================================

echo "...(8/11) SETTING UP SCRIPT LOOP"
echo " "
cat <<EOT >> snapshots/renbatchwnba.tcl

# II. PREPARE ORDER OF ITERATIION FOR STRUCTURE FILES
set mxx 2
set myy 3
for {set i 1} {AAi <= ${cn}} {incr i} { 

EOT

for ((f=1; f<=cn; f++)) 
 do
  echo "if {AAi=="$f"} {set cfile "AAfileid$f"}" >> snapshots/renbatchwnba.tcl
 done

# <=======================================================
# <==== IX. SETTING UP SCRIPT COMMANDS
# <=======================================================

echo "...(9/11) SETTING UP SCRIPT COMMANDS"
echo " "
cat <<EOT >> snapshots/renbatchwnba.tcl

# III. FORMAT REPRESENTATION AND DISPLAY
set j [expr {AAi - 1}]
mol new AAcfile.gro 
display depthcue off
mol modstyle "all" AAj VDW
mol modcolor "all" AAj Resname 
color Resname ${molec1} ${color1}
color Resname ${molec2} ${color2}
color Display Background white
scale by 1
scale to ${sclset}
display resize 1924 1061
display shadows on
display ambientocclusion on
axes location off
light 2 on
light 3 on

if {AAi==AAmxx} {
rotate y by 90
incr mxx 4
} 
if {AAi==AAmyy} {
rotate x by 90
incr myy 4
} 

# IV. RENDER IMAGES
puts "RENDERING FILE AAi of ${cn}: AAcfile.tga"
render Tachyon AAcfile "/usr/local/lib/vmd/tachyon_LINUXAMD64" -aasamples 12 %s -res 2000 2000 -format TARGA -o %s.tga 
puts "...DONE RENDERING"
puts "  "

# V. PREPARE FOR NEXT ITERATION
mol delete AAj
file delete AAcfile
}
exit
EOT

# <=======================================================
# <==== X. GENERATE SCRIPT FOR INDIVIDUAL RENDERING
# <=======================================================

echo "...(10/11) GENERATE SCRIPT FOR INDIVIDUAL RENDERING"
echo " "
cat <<EOT >> snapshots/renindivwnba.tcl
# WNBA CLUSTER SNAPSHOT INDIVIDUAL RENDERING SCRIPT 08/28/20
# CLUSTER SNAPSHOT IS RENDERED INTO TGA FILE WITH RES 2000 x 2000

# I. STRUCTURE FILE
set fileid 		# <==== insert structure file (without ext.) here
set cfile AAfileid

# II. FORMAT REPRESENTATION AND DISPLAY
mol new AAcfile.gro 
display depthcue off
mol modstyle "all" 0 VDW
mol modcolor "all" 0 Resname 
color Resname ${molec1} ${color1}
color Resname ${molec2} ${color2}
color Display Background white
scale by 1
scale to ${sclset}
display resize 1924 1061
display shadows on
display ambientocclusion on
axes location off
light 2 on
light 3 on

# III. RENDER IMAGES
puts "RENDERING FILE: indiv_AAcfile.tga"
render Tachyon indiv_AAcfile "/usr/local/lib/vmd/tachyon_LINUXAMD64" -aasamples 12 %s -res 2000 2000 -format TARGA -o %s.tga 
puts "  "
puts "...DONE RENDERING"
puts "  "

# IV. DELTE EXCESS FILE
file delete indiv_AAcfile
exit
EOT

vim -c ":%s/AA/$/g" -c ":wq" snapshots/renbatchwnba.tcl
vim -c ":%s/AA/$/g" -c ":wq" snapshots/renindivwnba.tcl


# <=======================================================
# <==== XI. REMOVE EXCESS FILES
# <=======================================================

echo "...(11/11) REMOVE EXCESS FILES"
echo " "
rm snapshots/temp2.txt
#rm *MER*

# =========================================================
# ------------------------- END ---------------------------
# =========================================================

echo "...DONE"
echo " "

echo "...FOR RENDERING, ENTER THE FOLLOWING COMMANDS"
echo " "
echo "   $ cd snapshots"
echo " "
echo "   $ vmd -dispdev none -e renbatchwnba.tcl"
echo " "
echo "   $ mogrify -antialias -density 300 -format png -quality 100 -colorspace RGB *.tga"
echo " "

duration=$SECONDS
echo " | TIME ELAPSED: $(($duration / 60)) MINUTE/S and $(($duration % 60)) SECOND/S |"
echo " "
