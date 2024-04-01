#!/bin/bash
# PART2 OF ANALYSIS SCRIPTS 
# DATE WRITTEN: AUGUST 27, 2020                                             
# DATE UPDATED: FEBRUARY 10, 2022
# PURPOSE: generate cluster reports

echo " "
echo " | COMMENCE: PART2 OF ANALYSIS SCRIPTS | GENERATE CLUSTER REPORTS | "
echo " "

# ========================================================
# ------------------    REMINDERS    ---------------------
# ========================================================

# 1. the sequence of analysis scripts are as follow:
#    [1] ./ANALYSIS_CLEAN.sh
#    [2] ./ANALYSIS_EXAMINE.sh   (*)
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
#    a. fort.10
#    b. cleaned trajectory file 
#
# 4. declare important variables in section I, for:
#    a. cleaned trajectory file to be read
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
# II.   SET SCRIPT FUNCTIONS
#
# III.  DETERMINE PRODUCTION PARAMETERS

# <=======================================================
# <==== I. DECLARE VARIABLES
# <=======================================================

echo "...(1/3) DECLARE VARIABLES"

clfile="system_mdd1-2_BUT_MET_TC500_clean.gro"

# <=======================================================
# <==== II. CONDITION CHECK
# <=======================================================

echo "...(2/3) CONDITION CHECK"

if compgen -G "fort.3*" > /dev/null; then
    echo " WARNING [1]: OUTPUT FILES PRESENT. TRANSFER BEFORE PROCEEDING" 
    echo " "
    exit 0
fi

# <=======================================================
# <==== III. PROBE SIMULATION DATA 
# <=======================================================

echo "...(3/3) PROBE SIMULATION DATA"

gfortran -o3 cluster_nm.f -o cluster || exit 0
./cluster < $clfile              

# =========================================================
# ------------------------- END ---------------------------
# =========================================================

echo "...DONE"
echo " "

duration=$SECONDS
echo " | TIME ELAPSED: $(($duration / 60)) MINUTE/S and $(($duration % 60)) SECOND/S |"
echo " "
