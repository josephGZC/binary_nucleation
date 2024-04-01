#!/bin/bash
# PART1 OF ANALYSIS SCRIPTS
# DATE WRITTEN: AUGUST 27, 2020
# DATE UPDATED: FEBRUARY 09, 2022
# PURPOSE: script to convert structure file to ca-readable file

echo " "
echo " | COMMENCE: PART1 OF ANALYSIS SCRIPTS | CONVERT STRUCTURE FILE | "
echo " "

# ========================================================
# ------------------    REMINDERS    ---------------------
# ========================================================

# 1. the sequence of analysis scripts are as follow:
#    [1] ./ANALYSIS_CLEAN.sh       (*)
#    [2] ./ANALYSIS_EXAMINE.sh   
#    [3] ./ANALYSIS_PLOT.sh
#    [4] ./ANALYSIS_SNAPSHOT.sh
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
#    a. trajectory file to be cleaned
#
# 4. declare important variables in section I, for:
#    a. trajectory file to be cleaned
#
# 5. if script is run in:
#    a. CSRC
#       nohup "./ANALYSIS_CLEAN.sh" > out.log &
#    b. ASTI
#       sbatch ANALYSIS_CLEAN.sub

# ========================================================
# -------------------    SECTIONS    ---------------------
# ========================================================

# I. DECLARE VARIABLES
#
# II. GENERATE CA-READABLE FILE

# <=======================================================
# <==== I. DECLARE VARIABLES
# <=======================================================

echo "...(1/2) DECLARE VARIABLES"
echo " "

basest="system_mdd1-5_SOL_BUT_TC500"
struct="${basest}.gro"                    #xtc trajectory to gro file
cleand="${basest}_clean.gro"

# <=======================================================
# <==== II. GENERATING CA-READABLE FILE
# <=======================================================

echo "...(2/2) GENERATING CA-READABLE FILE"   #cluster analysis readable
echo " "
SED_COM="s/^.\{5\}//;"                        #removes first 5 characters
SED_COM+="s/\(.\{3\}\).\{12\}/\1   /;"        #replaces 3rd to 12 character to white space
SED_COM+="p"                                  #print

sed '/[=]/d;/[A-Z]/!d' $struct | sed -n "$SED_COM" > $cleand    #all lines w/o capital letter are deleted
                                                                #sed -n means not to print. p tells sed to print what is matched
# =========================================================
# ------------------------- END ---------------------------
# =========================================================

echo "...DONE" 
echo " "

duration=$SECONDS
echo " | TIME ELAPSED: $(($duration / 60)) MINUTE/S and $(($duration % 60)) SECOND/S |"
echo " "
