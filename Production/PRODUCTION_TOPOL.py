#!/usr/bin/python3
# PART1 OF PRODUCTION SCRIPTS 09/06/20
# DATE WRITTEN: SEPTEMBER 5, 2020
# DATE UPDATED: OCTOBER 25, 2022
# PURPOSE: generate topol and mdp file

# <=======================================================
# <==== I.A. DECLARE VARIABLES (topol.top)
# <=======================================================

dec1        = "BUT"
dec2        = "MET"
topology1   = 'topol_unit_' + str(dec1) + '.top'
topology2   = 'topol_unit_' + str(dec2) + '.top'
topology3   = 'topol_system_' + str(dec1) + '_' + str(dec2) + '.top'
system      = 'unary alkanes'
molec1      = 'BUTA'
molec2      = 'MET'
molpop      = 10000

# <=======================================================
# <==== I.B. DECLARE VARIABLES (minim.mdp)
# <=======================================================

minimf      = 'minim.mdp'
emtol       = 10
nsteps1     = 12000

# <=======================================================
# <==== I.C. DECLARE VARIABLES (nvt.mdp)
# <=======================================================

nvtf        = 'nvt.mdp'
nsteps2     = 20000
temp1       = 700

# <=======================================================
# <==== I.D. DECLARE VARIABLES (md.mdp)
# <=======================================================

mdf         = 'md.mdp'
nsteps3     = 15000000
nst         = 500
temp2       = 230

# <=======================================================
# <==== II.A. GENERATE TOPOLOGY FILE - UNIT1
# <=======================================================

print("\nGenerating file: " + topology1)

with open(topology1, 'w+') as tp:
    prompt = 	'; WNBA Topology file '
    prompt +=   '\n\n#include "TraPPEwALKba.ff/forcefield.itp"'
    prompt +=   '\n#include "TraPPEwALKba.ff/tip4pew.itp"'
    prompt +=   '\n#include "TraPPEwALKba.ff/ammonia.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/ethane.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/propane.itp"'
    prompt +=	'\n#include "TraPPEwALKba.ff/butane.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/pentane.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/hexane.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/heptane.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/octane.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/nonane.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/decane.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/undecane.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/dodecane.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/butanol.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/methanol.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/octanol.itp"'
    prompt += 	'\n\n[ system ]'
    prompt += 	'\n' + system
    prompt += 	'\n\n[ molecules ]\n'
    tp.write(prompt)

with open(topology1, 'a') as tp:	
    prompt = molec1 + "   1\n" 			
    tp.write(prompt)

print("...Done!")

# <=======================================================
# <==== II.B. GENERATE TOPOLOGY FILE - UNIT2
# <=======================================================

print("\nGenerating file: " + topology2)

with open(topology2, 'w+') as tp:
    prompt = 	'; WNBA Topology file '
    prompt +=   '\n\n#include "TraPPEwALKba.ff/forcefield.itp"'
    prompt +=   '\n#include "TraPPEwALKba.ff/tip4pew.itp"'
    prompt +=   '\n#include "TraPPEwALKba.ff/ammonia.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/ethane.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/propane.itp"'
    prompt +=	'\n#include "TraPPEwALKba.ff/butane.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/pentane.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/hexane.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/heptane.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/octane.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/nonane.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/decane.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/undecane.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/dodecane.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/butanol.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/methanol.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/octanol.itp"'
    prompt += 	'\n\n[ system ]'
    prompt += 	'\n' + system
    prompt += 	'\n\n[ molecules ]\n'
    tp.write(prompt)

with open(topology2, 'a') as tp:	
    prompt = molec2 + "   1\n" 			
    tp.write(prompt)

print("...Done!")

# <=======================================================
# <==== II.C. GENERATE TOPOLOGY FILE - SYSTEM
# <=======================================================

print("\nGenerating file: " + topology3)

with open(topology3, 'w+') as tp:
    prompt = 	'; WNBA Topology file '
    prompt +=   '\n\n#include "TraPPEwALKba.ff/forcefield.itp"'
    prompt +=   '\n#include "TraPPEwALKba.ff/tip4pew.itp"'
    prompt +=   '\n#include "TraPPEwALKba.ff/ammonia.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/ethane.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/propane.itp"'
    prompt +=	'\n#include "TraPPEwALKba.ff/butane.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/pentane.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/hexane.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/heptane.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/octane.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/nonane.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/decane.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/undecane.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/dodecane.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/butanol.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/methanol.itp"'
    prompt += 	'\n#include "TraPPEwALKba.ff/octanol.itp"'
    prompt += 	'\n\n[ system ]'
    prompt += 	'\n' + system
    prompt += 	'\n\n[ molecules ]\n'
    tp.write(prompt)

count = 0
while count < int(molpop):
	count += 1
	with open(topology3, 'a') as tp:	
		if count % 2 == 0:		
			prompt = molec2 + "   1\n" 			
			tp.write(prompt)
		else:		
			prompt = molec1 + "   1\n" 			
			tp.write(prompt)
	if count == int(molpop):
		break

print("...Done!")

# <=======================================================
# <==== III. GENERATE MINIM.MDP FILE
# <=======================================================

print("\nGenerating file: " + minimf)

with open(minimf, 'w+') as mn:
    prompt =  'title					= minimization'
    prompt += '\nintegrator				= steep			; Algorithm (steep = steepest descent minimization)'
    prompt += '\nemtol					= ' + str(emtol) + '	; Stop minimization when the maximum force < 80.0 kJ/mol/nm'
    prompt += '\nemstep     				= 0.002     		; Energy step size, 2 fs'
    prompt += '\nnsteps					= ' + str(nsteps1) + '	; Maximum number of (minimization) steps to perform'
    prompt += '\n; Parameters describing how to find the neighbors of each atom and how to calculate the interactions'
    prompt += '\nnstlist	       			= 1		        ; Frequency to update the neighbor list and long range forces'
    prompt += '\ncutoff-scheme            		= Verlet'
    prompt += '\nns_type		       		= grid		    	; Method to determine neighbor list (simple, grid)'
    prompt += '\ncoulombtype				= cut-off		; Treatment of long range electrostatic interactions'
    prompt += '\npme_order	   			= 4			; cubic interpolation'
    prompt += '\n;Override fourierspacing, the best choice is powers of 2, 3, 5 and 7. Avoid large primes. For optimizing the relative load of the particle-particle interactions and the mesh part of PME, it is useful to know that the accuracy of the electrostatics remains nearly constant when the Coulomb cut-off (1.4 nm) and the PME grid spacing (50/32 = 1.56) are scaled by the same factor'
    prompt += '\nfourier-nx				= 32'
    prompt += '\nfourier-ny				= 32'
    prompt += '\nfourier-nz				= 32'
    prompt += '\nrcoulomb	        	        = 1.4		    	; Short-range electrostatic cut-off'
    prompt += '\nrvdw		         		= 1.4	    		; Short-range Van der Waals cut-off'
    prompt += '\npbc		  	    		= xyz			; Periodic Boundary Conditions (yes/no)'
    mn.write(prompt)

print("...Done!")

# <=======================================================
# <==== IV. GENERATE NVT.MDP FILE
# <=======================================================

print("\nGenerating file: " + nvtf)

with open(nvtf, 'w+') as nvt:
    prompt = 'title                                     = NVT equilibration' 
    prompt += '\nintegrator                             = md                            ; leap-frog integrator'
    prompt += '\nnsteps                                 = ' + str(nsteps2) + '          ; 2 * 5000 = 10000 fs = 10 ps'
    prompt += '\ndt                                     = 0.002                         ; 2 fs'
    prompt += '\n; Output control'
    prompt += '\nnstxout                                = 50                            ; save coordinates every 1.0 ps'
    prompt += '\nnstvout                                = 50                            ; save velocities every 1.0 ps'
    prompt += '\nnstenergy                              = 50                            ; save energies every 1.0 ps'
    prompt += '\nnstlog                                 = 50                            ; update log file every 1.0 ps'
    prompt += '\n; Bond parameters'
    prompt += '\ncontinuation                           = no                            ; first dynamics run'
    prompt += '\nconstraint_algorithm                   = lincs                         ; holonomic constraints' 
    prompt += '\nconstraints                            = all-bonds                     ; all bonds (even heavy atom-H bonds) constrained'
    prompt += '\nlincs_iter                             = 1		    		; accuracy of LINCS'
    prompt += '\nlincs_order                            = 4		    		; also related to accuracy'
    prompt += '\n; Neighborsearching'
    prompt += '\ncutoff-scheme                          = Verlet'
    prompt += '\nns_type		   	        = grid				; search neighboring grid cells'
    prompt += '\nnstlist		   	        = 10				; 20 fs, largely irrelevant with Verlet'
    prompt += '\nrcoulomb	   			= 1.4				; short-range electrostatic cutoff (in nm)'
    prompt += '\nrvdw		  			= 1.4				; short-range van der Waals cutoff (in nm)'
    prompt += '\n; Electrostatics'
    prompt += '\ncoulombtype				= PME				; Treatment of long range electrostatic interactions'
    prompt += '\npme_order	   			= 4		    		; cubic interpolation'
    prompt += '\n;fourierspacing			= 0.16				; grid spacing for FFT'
    prompt += '\n;Override fourierspacing, the best choice is powers of 2, 3, 5 and 7. Avoid large primes. For optimizing the relative load of the particle-particle interactions and the mesh part of PME, it is useful to know that the accuracy of the electrostatics remains nearly constant when the Coulomb cut-off (1.4 nm) and the PME grid spacing (50/32 = 1.56) are scaled by the same factor'
    prompt += '\nfourier-nx				= 32'
    prompt += '\nfourier-ny				= 32'
    prompt += '\nfourier-nz				= 32'
    prompt += '\n; Temperature coupling is on'
    prompt += '\ntcoupl					= V-rescale			; modified Berendsen thermostat'
    prompt += '\ntc-grps				= System			; two coupling groups - more accurate'
    prompt += '\ntau_t					= 0.1	        		; time constant, in ps'
    prompt += '\nref_t					= ' + str(temp1) + '	   	; reference temperature, one for each group, in K'
    prompt += '\n; Pressure coupling is off'
    prompt += '\npcoupl					= no 				; no pressure coupling in NVT'
    prompt += '\n; Periodic boundary conditions'
    prompt += '\npbc					= xyz				; 3-D PBC'
    prompt += '\n; Dispersion correction'
    prompt += '\nDispCorr				= EnerPres			; account for cut-off vdW scheme'
    prompt += '\n; Velocity generation'
    prompt += '\ngen_vel			    	= yes				; assign velocities from Maxwell distribution'
    prompt += '\ngen_temp				= ' + str(temp1) + '		; temperature for Maxwell distribution'
    prompt += '\ngen_seed				= -1				; generate a random seed'
    nvt.write(prompt)

print("...Done!")

# <=======================================================
# <==== V. GENERATE MD.MDP FILE
# <=======================================================

print("\nGenerating file: " + mdf)

with open(mdf, 'w+') as md:
    prompt =  'title                            = MD nucleation simulation'
    prompt += '\nintegrator                     = md                            ; leap-frog integrator'
    prompt += '\nnsteps                         = ' + str(nsteps3) + '           ; 2 * 1000000 = 2000000 fs (2 ns)'
    prompt += '\ndt                             = 0.002                         ; 2 fs'
    prompt += '\n; Output control'
    prompt += '\nnstxout                        = ' + str(nst) + '              ; save coordinates every 10.0 ps'
    prompt += '\nnstvout                        = ' + str(nst) + '              ; save velocities every 10.0 ps'
    prompt += '\nnstenergy                      = ' + str(nst) + '              ; save energies every 10.0 ps'
    prompt += '\nnstlog                         = ' + str(nst) + '              ; update log file every 10.0 ps'
    prompt += '\nnstxout-compressed             = ' + str(nst) + '              ; save compressed coordinates every 10.0 ps'
    prompt += '\n                                                               ; nstxout-compressed replaces nstxtcout'
    prompt += '\ncompressed-x-grps  		= System                        ; replaces xtc-grps'
    prompt += '\n; Bond parameters'
    prompt += '\ncontinuation	       		= yes                           ; Not restarting after NVT'
    prompt += '\nconstraint_algorithm           = lincs                         ; holonomic constraints'
    prompt += '\nconstraints	                = all-bonds                     ; all bonds (even heavy atom-H bonds) constrained'
    prompt += '\nlincs_iter	                = 1                             ; accuracy of LINCS'
    prompt += '\nlincs_order	                = 4                             ; also related to accuracy'
    prompt += '\n; Neighborsearching'
    prompt += '\ncutoff-scheme  	        = Verlet'
    prompt += '\nns_type		        = grid				; search neighboring grid cells'
    prompt += '\nnstlist		        = 10	   			; 20 fs, largely irrelevant with Verlet scheme'
    prompt += '\nrcoulomb	   	        = 1.4				; short-range electrostatic cutoff (in nm)'
    prompt += '\nrvdw		   	        = 1.4				; short-range van der Waals cutoff (in nm)'
    prompt += '\n; Electrostatics'
    prompt += '\ncoulombtype		        = PME				; Particle Mesh Ewald for long-range electrostatics'
    prompt += '\npme_order	    		= 4	   			; cubic interpolation'
    prompt += '\n;fourierspacing	        = 0.16				; grid spacing for FFT'
    prompt += '\nfourier-nx			= 32'
    prompt += '\nfourier-ny			= 32'
    prompt += '\nfourier-nz			= 32'
    prompt += '\n; Temperature coupling is on'
    prompt += '\ntcoupl		    	        = V-rescale	                ; modified Berendsen thermostat'
    prompt += '\ntc-grps			= System			; two coupling groups - more accurate'
    prompt += '\ntau_t			        = 0.1	                        ; time constant, in ps'
    prompt += '\nref_t			        = ' + str(temp2) + '             ; reference temperature, one for each group, in K'
    prompt += '\n; Pressure coupling is off in NVT'
    prompt += '\npcoupl			        = no'
    prompt += '\n;pcoupl		        = Parrinello-Rahman     	; Pressure coupling on in NPT'
    prompt += '\n;pcoupltype	     	        = isotropic	                ; uniform scaling of box vectors'
    prompt += '\n;tau_p		                = 2.0		                ; time constant, in ps'
    prompt += '\n;ref_p	    	                = 1.0		                ; reference pressure, in bar'
    prompt += '\n;compressibility  	        = 4.5e-5	                ; isothermal compressibility of water, bar^-1'
    prompt += '\n; Periodic boundary conditions'
    prompt += '\npbc			        = xyz	                        ; 3-D PBC'
    prompt += '\n; Dispersion correction'
    prompt += '\nDispCorr		        = EnerPres	                ; account for cut-off vdW scheme'
    prompt += '\n; Velocity generation'
    prompt += '\ngen_vel		        = no		                ; Velocity generation is off'
    md.write(prompt)

print("...Done!")

# ========================================================
# ------------------      END      -----------------------                  
# ========================================================
