### CONVERSION PARAMETERS (Xs slope + Xi intercept) ###
Ts="140.5678"    ### convert Temperature to K
Ti="185.4676"
Ps="469.459"     ### convert Pressure to MPa
Pi="-211.15552"
rhos="1527.3"    ### convert Density to kg/m^3
rhoi="-23.1017"
CPs="4.455"    ### convert Specific heat to J / g mol
CPi="3.56843"  
##
Hs="11.272132"     # Convert Enthalpy to kJ/mol (0.018*Ts*CPs)
Hi="9.028924"      #                            (0.018*Ts*CPi)
###
###
L="32"
for size in $L
do
for file in $(ls L$size\_P*_density.dat)
do
awk -v Ts=$Ts -v Ti=$Ti -v rhos=$rhos -v rhoi=$rhoi '{print $1*Ts+Ti, $2*rhos+rhoi, $3*rhos}' $file > $file.SI
done
for file in $(ls L$size\_P*_enthalpy.dat)
do
awk -v Ts=$Ts -v Ti=$Ti -v Hs=$Hs -v Hi=$Hi '{print $1*Ts+Ti, $2*Hs+Hi*$1, $3*Hi}' $file > $file.SI
done
for file in $(ls L$size\_P*_nhb.dat L32_P*_ncoop.dat)
do
awk -v Ts=$Ts -v Ti=$Ti '{print $1*Ts+Ti, $2, $3}' $file > $file.SI
done
for file in $(ls L$size\_P*_Cpfluc.dat)
do
awk -v Ts=$Ts -v Ti=$Ti -v CPs=$CPs -v CPi=$CPi '{print $1*Ts+Ti, $2*CPs+CPi, $3*CPs}' $file > $file.SI
done
### Transform isochores
#for file in $(ls L$size\_isocora_rho*.dat)
#do
#awk -v Ts=$Ts -v Ti=$Ti -v Ps=$Ps -v Pi=$Pi '{print $1*Ts+Ti, $2*Ps+Pi}' $file > $file.SI
#done
### Transform loci of extrema (phase diagram)
#for file in $(ls L$size\_J0.5_Js0.08_vHB0.6_*.dat)
#do
#awk -v Ts=$Ts -v Ti=$Ti -v Ps=$Ps -v Pi=$Pi '{print $1*Ts+Ti, $2*Ps+Pi, $3*Ts}' $file > $file.SI
#done
done
##############
##############  Alpha_P calculated as 1/V (\partial V/ \partial T)_P
##############
L="32"
P="-0.8 -0.75 -0.7 -0.6 -0.5 -0.4 -0.38 -0.35 -0.3 -0.25 -0.22 -0.2 -0.15 -0.13 -0.1 -0.05 0.0 0.05 0.1 0.13 0.15 0.2 0.24 0.25 0.3 0.343 0.35 0.38 0.4 0.45 0.47 0.5 0.515 0.53 0.55 0.6 0.63 0.66 0.7 0.73 0.77 0.8 0.85 0.88 0.9 0.95 1"
for size in $L
do
for pres in $P
do
####
    file=L$size\_P$pres\_density.dat.SI
    if [ -f $file ]; then 
      awk '{print $16}' L$size\_P$pres.dat > ndata.aux
      awk '{print $1, 1.0/$2 , $3*$3/(3*$2*$2*$2*$2)}' $file > dens.aux  #fluc. specific volume
      paste dens.aux ndata.aux > AlphaP_input.aux
      awk '{print $1, $2 , $3*$4, $4}' AlphaP_input.aux > AlphaP_input
      rm ndata.aux dens.aux AlphaP_input.aux
      ./AlphaPderivate   #Calculates Alpha_P as derivative
      mv AlphaP_derivate.out L$size\_P$pres\_AlphaPderivate.dat.SI
      awk '$2 != "-nan" && $2 != "inf" && $2 != "-inf" && $3 != "inf" {print $1, $2*1000, $3*1000}' L$size\_P$pres\_AlphaPderivate.dat.SI > L$size\_P$pres\_AlphaPderivate.dat.SI.aux  ## converts Alpha_P to 10^3 K^-1
      mv L$size\_P$pres\_AlphaPderivate.dat.SI.aux L$size\_P$pres\_AlphaPderivate.dat.SI
      rm AlphaP_input
    fi
##############
##############  Energy in SI units calculated as E = H - P/rho
##############
    file1=L$size\_P$pres\_enthalpy.dat.SI   # kJ/mol
    file2=L$size\_P$pres\_density.dat.SI    # kg/m^3
    if [ -f $file1 ]; then
      if [ -f $file2 ]; then
        paste $file1 $file2 > energy.aux
        awk -v P=$pres -v Ps=$Ps -v Pi=$Pi '{print $1, $2-18*(P*Ps+Pi)/$5, $3-18*(P*Ps+Pi)*$6/($5*$5)}' energy.aux >  L$size\_P$pres\_energy.dat.SI
        rm energy.aux
      fi
    fi  
done
done
##############
##############  K_T calculated as -1/V (\partial V/ \partial P)_T
##############
for s1 in $L
do
if [ -f KT_list_isobars ]; then
  rm KT_list_isobars
fi
if [ -f KT_list_isotherms ]; then
  rm KT_list_isotherms
fi
#### From L32_average_data.dat ; generate a file containing T, P, rho, error(rho) in SI units
awk -v Ts=$Ts -v Ti=$Ti -v Ps=$Ps -v Pi=$Pi -v rhos=$rhos -v rhoi=$rhoi '{print Ts*$3+Ti, Ps*$2+Pi, (rhos*($1*$1*$1)/$8+rhoi), (rhos*($1*$1*$1)*sqrt(3*$9/$16)/($8*$8))}' L$s1\_average_data.dat > KT_input
for p1 in $P
do
  echo $p1 >> KT_list_isobars
done
awk -v Ps=$Ps -v Pi=$Pi '{print Ps*$1+Pi}' KT_list_isobars > KT_list_isobars.aux
mv KT_list_isobars.aux KT_list_isobars  #Generate a list of isobars in MPa
#
# Since we need to derivate along isotherms, I generate also a list of temperatures, in K
TEMP="0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2"
for t1 in $TEMP
do
  echo $t1 >> KT_list_isotherms
done
awk -v Ts=$Ts -v Ti=$Ti '{print Ts*$1+Ti}' KT_list_isotherms > KT_list_isotherms.aux
mv KT_list_isotherms.aux KT_list_isotherms
./KTderivate #Calculates KT as derivative along isotherms and stores the result along isobars
rm KT_input KT_list_isobars KT_list_isotherms
rm T*_density_isoth.dat  #comment to keep density along isotherms
for file in $(ls P*KTderivate.dat)
do
 awk  '{print $1, $2*1000, $3*1000, $4}' $file > L$s1\_$file.SI  # convert to 10^3 MPa^-1 = GPa^-1
  sort -g L$s1\_$file.SI > L$s1\_$file.SI.aux
  mv L$s1\_$file.SI.aux L$s1\_$file.SI
  rm $file 
done
done
