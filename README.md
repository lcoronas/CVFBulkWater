# FSBulkWater
######################
####
#### A transferable molecular model for accurate thermodynamic studies of water in large-scale systems
####
#### Luis E. Coronas, Oriol Vilanova, and Giancarlo Franzese
####
#### Universitat de Barcelona
####
######################

This repository contains all the data to reproduce FS data in the main text and supplementary materials.

###########################
##### BULK FS MODEL SOURCE CODE
###########################

SOURCE folder contains the code and a sample of input/output files

A) gpu_water.cu : source code

   how to compile: nvcc -arch=sm_86 --ptxas-options=--verbose --use_fast_math -O2 -o gpu_water3D gpu_water3D.cu mersenne_inline.cu -lm
	where -arch=sm_86 should change according to "CUDA Capability Major/Minor version number"
	
B) input files:

	1) input_data contains
		Pressure                    0.45    	# 0.101 MPa
		Temperature                 0.815    	# 300 K
		Seed                        44648604	# Random number generator seed
		Metropolis_Steps            1		# Number of Metropolis updates per MC Step
		Cluster_Steps               0           # Number of Swendsen-Wang updates per MC Step
		Equilibrium_Steps           10000	# Equilibration steps (not recorded for statictics)
		Sampling_Steps              100000	# Number of steps (for production)
		Sampling_Interval           10		# Sample evrey x steps
		R_cutoff		    6.0		# Cutoff distance van der Waals
		Flag_config                 1		# Start from: 0 random config; 1 input config file; 2 fully ordered config
		Flag_distances              0		# Set to 0. If 1, uses pre-calculated distances in distances file
		Flag_correlation            0		# 0: Do not calculate autocorrelation function. 1 to enable
		Flag_chessboard             1		# Set to 1. Enabels checherboard algorithm for eta variables [http://hdl.handle.net/2445/197881]
		Device_ID                   0		# GPU ID. If multiple GPUs are available, enables choosing in which to run
		
	2) input_config
		For every molecule: state of the 6 sigma variables + state of the 6 eta variables
		
C) Output files

	1) water3D.log : logfile
	2) data.out : stores timeseries of thermodynamic observables
	    $1 time; $2 energy E; $3 enthalpy H; $4 volume V; $5 number of hydrogen bonds N_hb, $6 number of cooperative bonds N_sigma/N_sigma^Max, $7 Radius of the cells 
	3) cluster.out : stores timeseries of cluster analysis. Only valid if using cluster MC, not Metropolis.
	    $1 time; $2 largest_cluster_size, $3 second_largest_cluster_size, $4 num_percolating_clusters, $5 average cluster_size , $6 sigma^2(cluster_size), $7 num_clusters
	4) config_out : final configuration of the system (to use in sequential annealing or to resume the simulation)

############################
##### DATA AND ANALYSIS SCRITPS
###########################

FS_WATER_32x32x32_DATA folder contains rough data and analysis script.

A) WATER32_DATA.tar.gz : a tarred file containing L32_T*_P* folders containing FS simulation files [AVAILABLE ON REQUEST]

B) Analysis scripts and code

To analyze the data:

   1st) untar WATER32_DATA.tar.gz
   2nd) execute "bash script_analisi.sh"
   
   You may need to compile KTderivate.c, AlphaPderivate.c, and CPderivate.c
   Execute: gcc -o KTderivate KTderivate.c -lm
            gcc -o AlphaPderivate AlphaPderivate.c -lm
            gcc -o CPderivate CPderivate.c -lm
            
The script generates L32_P*_*.dat files [INTERNAL UNITS] and L32_P*_*.dat.SI files [INTERNATIONAL SYSTEM UNITS]
A description of the files used in the paper is in the following:

###########################
##### PROCESSED_DATA
###########################

PROCESSED_DATA folder contains already-processed data

A) FS_isobars.tar.gz contains 

	1) L32_P*_density.dat, L32_P*_enthalpy.dat : Density and enthalpy along isobar P [in internal units]

	$1 Temperature   [internal units]
	$2 Density / Enthalpy [internal units]
	$3 Error in density / enthalpy [internal units]

	2) L32_P*_density.dat.SI, L32_P*_AlphaPderivate.dat.SI, L32_Cpfluc.dat.SI : Density, Alpha_P and C_P along isobar P [in internal units]

	$1 Temperature  [K]
	$2 Density / Alpha_P / C_P [kg/m^3 ; K^-1 ; J g^-1 K^-1] 
	$3 Error in Density / Alpha_P / C_P [kg/m^3 ; K^-1 ; J g^-1 K^-1] 

	3) L32_P*_KTderivate.dat.SI : Isothermal compressibility along isobar P [MPa]
	$1 Temperature [K]
	$2 Isothermal compressibility [GPa^-1]
	$3 Error in Isothermal compressibility [GPa^-1]
	
B) FS_isochores.tar.gz contains 

	1) L32_isocora_rho*.dat [internal units]
	$1 Temperature [internal units]
	$2 Pressure [internal units]
	
	2) L32_isocora_rho*.dat.SI [internal units]
	$1 Temperature [K]
	$2 Pressure [MPa]
	
C) FS_phase_diagram.tar.gz contains

	1) L32_J0.5_Js0.08_vHB0.6_*.dat / L32_J0.5_Js0.08_vHB0.6_*.dat.SI : TMD, TminD, LGSpinodal, CPmax, CPmin, KTmin, KTmax, AlphaP=0 projections into Temperature-Pressure phase diagram
	$1 Temperature [internal units / K]
	$2 Pressure [internal units / MPa]
	$3 Temperautre error [internal units / K]
	
	
	2) L32_J0.5_Js0.08_vHB0.6_*_Density_T.dat /  L32_J0.5_Js0.08_vHB0.6_*_Density_T.dat.SI : TMD, TminD, LG Spindoal projection into Temperature-Density phase diagram
	$1 Temperature [internal units / K]
	$2 Density [internal units / kg m^-3]
	$3 Temperautre error [internal units / K]

 D) Pressure_conversion.dat contains list of isobars in internernal units and SI units
	
	
