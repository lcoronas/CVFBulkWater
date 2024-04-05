#!/bin/bash

#This scripts performs thermodynamic averages and error calculations on L32_T*_P* folders of Bulk FS water simulations

# It runs over all the folders and calculate the average values of the observable [INTERNAL UNITS].
# It generates a file "L32_average_data.dat" with 16 columns representing respectively:
# 1) Lattice Size L
# 2) Pressure
# 3) Temperature
# 4) Average Energy
# 5) Energy Fluctuation
# 6) Average Enthalpy
# 7) Enthalpy Fluctuation
# 8) Average Volume
# 9) Volume Fluctuation
# 10) Average Number of HB
# 11) Number of HB Fluctuation
# 12) Average Number of ( Cooperative Bonds/(15*N) )
# 13) Number of ( Cooperative Bonds/(15*N) ) Fluctuation
# 14) Average Cell Radius
# 15) Cell Radius Fluctuation
# 16) Number of Data

if [ -f L32_average_data.dat ]; 
then
  rm L32_average_data.dat
fi

for d in $(ls -d L32_T*_P*); 
do
  L=$(echo $d | grep -o -E '[0-9.-]+' | head -1)              
  T=$(echo $d | grep -o -E '[0-9.-]+' | sed '1d' | head -1)
  P=$(echo $d | grep -o -E '[0-9.-]+' | sed '1,2d' | head -1)
  echo $L $T $P
  awk -v L=$L -v T=$T -v P=$P 'BEGIN { OFS = "\t"}{a+=$2;a2+=$2*$2;b+=$3;b2+=$3*$3;c+=$4;c2+=$4*$4;d+=$5;d2+=$5*$5;e+=$6;e2+=$6*$6;f+=$7;f2+=$7*$7}END{print L,P,T,a/NR,a2/NR-(a/NR)*(a/NR),b/NR,b2/NR-(b/NR)*(b/NR),c/NR,c2/NR-(c/NR)*(c/NR),d/NR,d2/NR-(d/NR)*(d/NR),e/NR,e2/NR-(e/NR)*(e/NR),f/NR,f2/NR-(f/NR)*(f/NR),NR}' $d/data.out >> L32_average_data.dat
done
bash filter_isobars.sh  #Filters L32_average_data.dat into single-file per isobar
bash calculadora.sh     #Calculates thermodynamic averages over isobars [INTERNAL UNITS]
bash conversor.sh	#Converts thermodynamic averages into SI units and calculates Alpha_P and K_T from volume derivatives
