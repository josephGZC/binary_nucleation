#!/bin/bash
# PART3 OF ANALYSIS SCRIPTS 
# DATE WRITTEN: AUGUST 27, 2020                                             
# DATE UPDATED: FEBRUARY 10, 2022
# PURPOSE: generate graphs
SECONDS=0

echo " "
echo " | COMMENCE: PART3 OF ANALYSIS SCRIPTS | GENERATE GRAPHS | "
echo " "

# ========================================================
# ------------------    REMINDERS    ---------------------
# ========================================================

# 1. the sequence of analysis scripts are as follow:
#    [1] ./ANALYSIS_CLEAN.sh
#    [2] ./ANALYSIS_EXAMINE.sh   
#    [3] ./ANALYSIS_PLOT.sh       (*)
#    [4] ./ANALYSIS_SNAPSHOT.sh
#
# 2. further breakdown of the scripts are as follow:
#    [1] ./ANALYSIS_CLEAN.sh (insert structure file without ext.)
#
#    [2] ./ANALYSIS_EXAMINE.sh (insert clean structure file with ext.)
#         (A) cluster_nm.f
#
#    [3] ./ANALYSIS_PLOT.sh (insert clean structure file with ext.)
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
#    [4] ./ANALYSIS_SNAPSHOT.sh  (adjust scale variable sclset)
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
#    c. fort output files
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
# II.   GENERATE GRAPH-COMPATIBLE FILES
#
# III.  GENERATE RADIAL DENSITY PROFILE #1A
#
# IV.   GENERATE RADIAL DENSITY PROFILE #1B
#
# V.    GENERATE RADIAL DENSITY PROFILE #2A
#
# VI.   GENERATE RADIAL DENSITY PROFILE #2B
#
# VII.  GENERATE RADIAL DENSITY PROFILE #3A
#
# VIII. GENERATE RADIAL DENSITY PROFILE #3B
#
# IX.   GENERATE RADIAL DISTRIBUTION FUNCTION #1A
#
# X.    GENERATE RADIAL DISTRIBUTION FUNCTION #1B
#
# XI.   GENERATE FREQUENCY GRAPH
#
# XII.  GENERATE MOLE FRACTION GRAPH
#
# XIII. GENERATE MASS GRAPH
#
# XIV.  GENERATE MASS FRACTION (/MAX) GRAPH
#
# XV.   GENERATE MASS FRACTION (/MIN) GRAPH
#
# XVI.  GENERATE RADIUS GRAPH
#
# XVII. CREATE GRAPHS FOLDER

# <=======================================================
# <==== I. DECLARE VARIABLES 
# <=======================================================

echo "...(1/17) DECLARE VARIABLES"

clfile="system_mdd1-2_BUT_MET_TC500_clean.gro"

# <=======================================================
# <==== II. GENERATE GRAPH-COMPATIBLE FILES 
# <=======================================================

echo "...(2/17) GENERATE GRAPH-COMPATIBLE FILES"
echo " "

gfortran -o3 graphing.f -o graphing
./graphing < $clfile

sed -i 's/ \{3,\}//g' PLOTRDP1A.gnu         #removes consecutive spaces that are more than three
sed -i 's/ \{3,\}//g' PLOTRDP1B.gnu
sed -i 's/ \{3,\}//g' PLOTRDP2A.gnu
sed -i 's/ \{3,\}//g' PLOTRDP2B.gnu
sed -i 's/ \{3,\}//g' PLOTRDP3A.gnu
sed -i 's/ \{3,\}//g' PLOTRDP3B.gnu
sed -i 's/ \{3,\}//g' PLOTRDF1A.gnu
sed -i 's/ \{3,\}//g' PLOTRDF1B.gnu
sed -i 's/ \{3,\}//g' PLOTmolfraction.gnu

# <=======================================================
# <==== III. GENERATE RADIAL DENSITY PROFILE #1A
# <=======================================================

if [[ -f "PLOTRDP1A.gnu" ]]; then
echo "...(3/17) GENERATE RADIAL DENSITY PROFILE #1A"
echo " "
./PLOTRDP1A.gnu

while read pl; do
if [[ $pl =~ output ]]; then
nn="${pl:12:-5}.txt" 
elif [[ $pl =~ "plot 'fort" ]]; then
nd=${pl:6:9}
cp $nd $nn
fi
done < PLOTRDP1A.gnu

else
echo "...(3/17) SKIPPING CREATION OF RADIAL DENSITY PROFILE #1A"
echo " "
fi

# <=======================================================
# <==== IV. GENERATE RADIAL DENSITY PROFILE #1B
# <=======================================================

if [[ -f "PLOTRDP1B.gnu" ]]; then
echo "...(4/17) GENERATE RADIAL DENSITY PROFILE #1B"
echo " "
./PLOTRDP1B.gnu

while read pl; do
if [[ $pl =~ output ]]; then
nn="${pl:12:-5}.txt" 
elif [[ $pl =~ "plot 'fort" ]]; then
nd=${pl:6:9}
cp $nd $nn
fi
done < PLOTRDP1B.gnu

else
echo "...(4/17) SKIPPING CREATION OF RADIAL DENSITY PROFILE #1B"
echo " "
fi

# <=======================================================
# <==== V. GENERATE RADIAL DENSITY PROFILE #2A
# <=======================================================

if [[ -f "PLOTRDP2A.gnu" ]]; then
echo "...(5/17) GENERATE RADIAL DENSITY PROFILE #2A"
echo " "
./PLOTRDP2A.gnu 

while read pl; do
if [[ $pl =~ output ]]; then
nn="${pl:12:-5}.txt" 
elif [[ $pl =~ "plot 'fort" ]]; then
nd=${pl:6:9}
cp $nd $nn
fi
done < PLOTRDP2A.gnu

else
echo "...(5/17) SKIPPING CREATION OF RADIAL DENSITY PROFILE #2A"
echo " "
fi  

# <=======================================================
# <==== VI. GENERATE RADIAL DENSITY PROFILE #2B
# <=======================================================

if [[ -f "PLOTRDP2B.gnu" ]]; then
echo "...(6/17) GENERATE RADIAL DENSITY PROFILE #2B"
echo " "
./PLOTRDP2B.gnu 

while read pl; do
if [[ $pl =~ output ]]; then
nn="${pl:12:-5}.txt" 
elif [[ $pl =~ "plot 'fort" ]]; then
nd=${pl:6:9}
cp $nd $nn
fi
done < PLOTRDP2B.gnu

else 
echo "...(6/17) SKIPPING CREATION OF RADIAL DENSITY PROFILE #2B"
echo " "
fi

# <=======================================================
# <==== VII. GENERATE RADIAL DENSITY PROFILE #3A
# <=======================================================

if [[ -f "PLOTRDP3A.gnu" ]]; then
echo "...(7/17) GENERATE RADIAL DENSITY PROFILE #3A"
echo " "
./PLOTRDP3A.gnu 

while read pl; do
if [[ $pl =~ output ]]; then
nn="${pl:12:-5}.txt" 
elif [[ $pl =~ "plot 'fort" ]]; then
nd=${pl:6:9}
cp $nd $nn
fi
done < PLOTRDP3A.gnu

else
echo "...(7/17) SKIPPING CREATION OF RADIAL DENSITY PROFILE #3A"
echo " "
fi

# <=======================================================
# <==== VIII. GENERATE RADIAL DENSITY PROFILE #3B
# <=======================================================

if [[ -f "PLOTRDP3B.gnu" ]]; then
echo "...(8/17) GENERATE RADIAL DENSITY PROFILE #3B"
echo " "
./PLOTRDP3B.gnu 

while read pl; do
if [[ $pl =~ output ]]; then
nn="${pl:12:-5}.txt" 
elif [[ $pl =~ "plot 'fort" ]]; then
nd=${pl:6:9}
cp $nd $nn
fi
done < PLOTRDP3B.gnu

else
echo "...(8/17) SKIPPING CREATION OF RADIAL DENSITY PROFILE #3B"
echo " "
fi

# <=======================================================
# <==== IX. GENERATE RADIAL DISTRIBUTION FUNCTION #1A
# <=======================================================

if [[ -f "PLOTRDF1A.gnu" ]]; then
echo "...(9/17) GENERATE RADIAL DISTRIBUTION FUNCTION #1A"
echo " "
./PLOTRDF1A.gnu

while read pl; do
if [[ $pl =~ output ]]; then
nn="${pl:12:-5}.txt" 
elif [[ $pl =~ "plot 'fort" ]]; then
nd=${pl:6:9}
cp $nd $nn
fi
done < PLOTRDF1A.gnu

else
echo "...(9/17) SKIPPING CREATION OF RADIAL DISTRIBUTION FUNCTION #1A"
echo " "
fi

# <=======================================================
# <==== X. GENERATE RADIAL DISTRIBUTION FUNCTION #1B
# <=======================================================

if [[ -f "PLOTRDF1B.gnu" ]]; then
echo "...(10/17) GENERATE RADIAL DISTRIBUTION FUNCTION #1B"
echo " "
./PLOTRDF1B.gnu

while read pl; do
if [[ $pl =~ output ]]; then
nn="${pl:12:-5}.txt" 
elif [[ $pl =~ "plot 'fort" ]]; then
nd=${pl:6:9}
cp $nd $nn
fi
done < PLOTRDF1B.gnu

else
echo "...(10/17) SKIPPING CREATION OF RADIAL DISTRIBUTION FUNCTION #1B"
echo " "
fi

# <=======================================================
# <==== XI. GENERATE FREQUENCY GRAPH
# <=======================================================

if [[ -f "PLOTfrequency.gnu" ]]; then
echo "...(11/17) GENERATE FREQUENCY GRAPH"
echo " "
./PLOTfrequency.gnu     

while read pl; do
if [[ $pl =~ output ]]; then
nn="${pl:12:-5}.txt" 
elif [[ $pl =~ "plot 'fort" ]]; then
nd=${pl:7:9}
cp $nd $nn
fi
done < PLOTfrequency.gnu

else
echo "...(11/17) SKIPPING CREATION OF FREQUENCY GRAPH"
echo " "
fi

# <=======================================================
# <==== XII. GENERATE MOLE FRACTION GRAPH
# <=======================================================

if [[ -f "PLOTmolfraction.gnu" ]]; then
echo "...(12/17) GENERATE MOLE FRACTION GRAPH"
echo " "
./PLOTmolfraction.gnu 

while read pl; do
if [[ $pl =~ output ]]; then
nn="${pl:12:-5}.txt" 
elif [[ $pl =~ "plot 'fort" ]]; then
nd=${pl:7:9}
cp $nd $nn
fi
done < PLOTmolfraction.gnu

else
echo "...(12/17) SKIPPING CREATION OF MOLE FRACTION GRAPH"
echo " "
fi

# <=======================================================
# <==== XIII. GENERATE MASS GRAPH
# <=======================================================

if [[ -f "PLOTmass.gnu" ]]; then
echo "...(13/17) GENERATE MASS GRAPH"
echo " "
./PLOTmass.gnu  

while read pl; do
if [[ $pl =~ output ]]; then
nn="${pl:12:-5}.txt" 
elif [[ $pl =~ "plot 'fort" ]]; then
nd=${pl:7:9}
cp $nd $nn
fi
done < PLOTmass.gnu

else
echo "...(13/17) SKIPPING CREATION OF MASS GRAPH"
echo " "
fi

# <=======================================================
# <==== XIV. GENERATE MASS FRACTION (/MAX) GRAPH
# <=======================================================

if [[ -f "PLOTmassfracmax.gnu" ]]; then
echo "...(14/17) GENERATE MASS FRACTION (/MAX) GRAPH"
echo " "
./PLOTmassfracmax.gnu

while read pl; do
if [[ $pl =~ output ]]; then
nn="${pl:12:-5}.txt" 
elif [[ $pl =~ "plot 'fort" ]]; then
nd=${pl:7:9}
cp $nd $nn
fi
done < PLOTmassfracmax.gnu

else 
echo "...(14/17) SKIPPING CREATION OF MASS FRACTION (/MAX) GRAPH"
echo " "
fi

# <=======================================================
# <==== XV. GENERATE MASS FRACTION (/MIN) GRAPH
# <=======================================================

if [[ -f "PLOTmassfracmin.gnu" ]]; then
echo "...(15/17) GENERATE MASS FRACTION (/MIN) GRAPH"
echo " "
./PLOTmassfracmin.gnu 

while read pl; do
if [[ $pl =~ output ]]; then
nn="${pl:12:-5}.txt" 
elif [[ $pl =~ "plot 'fort" ]]; then
nd=${pl:7:9}
cp $nd $nn
fi
done < PLOTmassfracmin.gnu

else
echo "...(15/17) SKIPPING CREATION OF MASS FRACTION (/MIN) GRAPH"
echo " "
fi

# <=======================================================
# <==== XVI. GENERATE RADIUS GRAPH
# <=======================================================

if [[ -f "PLOTradius.gnu" ]]; then
echo "...(16/17) GENERATE RADIUS GRAPH"
echo " "
./PLOTradius.gnu

while read pl; do
if [[ $pl =~ output ]]; then
nn="${pl:12:-5}.txt" 
elif [[ $pl =~ "plot 'fort" ]]; then
nd=${pl:7:9}
cp $nd $nn
fi
done < PLOTradius.gnu

else
echo "...(16/17) SKIPPING CREATION OF RADIUS GRAPH"
echo " "
fi

# <=======================================================
# <==== XVII. CREATE GRAPHS FOLDER
# <=======================================================

echo "...(17/17) EXPORTING GRAPHS TO PNG FORMAT"
echo " "

mogrify -antialias -density 300 -format png -quality 100 -colorspace RGB *.eps

[[ ! -d "graphs" ]] &&
mkdir graphs &&
mkdir graphs/eps &&
mkdir graphs/png &&
mkdir graphs/gnu &&
mkdir graphs/txt

mv *.eps graphs/eps 
mv *.png graphs/png
mv *.gnu graphs/gnu
mv *.txt graphs/txt
    
# =========================================================
# ------------------------- END ---------------------------
# =========================================================

echo "...DONE"
echo " "

duration=$SECONDS
echo " | TIME ELAPSED: $(($duration / 60)) MINUTE/S and $(($duration % 60)) SECOND/S |"
echo " "
