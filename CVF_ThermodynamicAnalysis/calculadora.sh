# This script generates files L32_P*_f.dat where f is a thermodynamic observable
# f: energy, enthalpy, volume, nhb, ncoop.
# Each file contains:
#     $1 Temperature
#     $2 Average of f: <f> = <F>/N
#     $3 Error of the average: error(<f>) = sqrt[ (1+2*tau/deltaT)*sigma^2(F) /n_measurements ]/N ~ sqrt[ 3*sigma^2(F) /n_measurements ]/N ; where tau is the correlation time and deltaT the sampling time interval
# Note: density = 1/V ; error = error(V)/V^2
# Note: ncoop is already normalized to Ncoop^Max
#
# L32_P*_Cpfluc.dat
#     $1 Temperature
#     $2 Cp = sigma^2(H)/(T^2 N)

L="32"
P="-0.8 -0.75 -0.7 -0.6 -0.5 -0.4 -0.38 -0.35 -0.3 -0.25 -0.22 -0.2 -0.15 -0.13 -0.1 -0.05 0.0 0.05 0.1 0.13 0.15 0.2 0.24 0.25 0.3 0.343 0.35 0.38 0.4 0.45 0.47 0.5 0.515 0.53 0.55 0.6 0.63 0.66 0.7 0.73 0.77 0.8 0.85 0.88 0.9 0.95 1"
######
for size in $L
do
for pres in $P
do
  file=L$size\_P$pres.dat  
  awk '{print $3, $4/($1*$1*$1), sqrt(3*$5/$16)/($1*$1*$1)}' $file > L$size\_P$pres\_energy.dat
     sort -g L$size\_P$pres\_energy.dat > L$size\_P$pres\_energy.aux
     mv L$size\_P$pres\_energy.aux L$size\_P$pres\_energy.dat
  awk '{print $3, $6/($1*$1*$1), sqrt(3*$7/$16)/($1*$1*$1)}' $file > L$size\_P$pres\_enthalpy.dat
  awk '{print $3, $8/($1*$1*$1), sqrt(3*$9/$16)/($1*$1*$1)}' $file > L$size\_P$pres\_volume.dat
  awk '{print $3, ($1*$1*$1)/$8, ($1*$1*$1)*sqrt(3*$9/$16)/($8*$8)}' $file > L$size\_P$pres\_density.dat
     sort -g L$size\_P$pres\_density.dat > L$size\_P$pres\_density.aux
     mv L$size\_P$pres\_density.aux L$size\_P$pres\_density.dat
  awk '{print $3, $10/(2*$1*$1*$1), sqrt(3*$11/$16)/(2*$1*$1*$1)}' $file > L$size\_P$pres\_nhb.dat
     sort -g L$size\_P$pres\_nhb.dat > L$size\_P$pres\_nhb.aux
     mv L$size\_P$pres\_nhb.aux L$size\_P$pres\_nhb.dat
  awk '{print $3, $12, sqrt(3*$13/$16)}' $file > L$size\_P$pres\_ncoop.dat
     sort -g L$size\_P$pres\_ncoop.dat > L$size\_P$pres\_ncoop.aux
     mv L$size\_P$pres\_ncoop.aux L$size\_P$pres\_ncoop.dat
  awk '{print $3, $14, sqrt($15/$16)}' $file > L$size\_P$pres\_rcell.dat
     sort -g L$size\_P$pres\_rcell.dat > L$size\_P$pres\_rcell.aux
     mv L$size\_P$pres\_rcell.aux L$size\_P$pres\_rcell.dat
  awk '{print $3, $7/($3*$3*$1*$1*$1)}' $file > L$size\_P$pres\_Cpfluc.dat
     sort -g L$size\_P$pres\_Cpfluc.dat > L$size\_P$pres\_Cpfluc.aux
     mv L$size\_P$pres\_Cpfluc.aux L$size\_P$pres\_Cpfluc.dat
  awk '{print $3, $9/($3*$8)}' $file > L$size\_P$pres\_KT.dat
  awk '{print $3, $8, $9, $16}' $file > AlphaP_input
  ./AlphaPderivate
  mv AlphaP_derivate.out L$size\_P$pres\_AlphaP.dat
  rm AlphaP_input
 # awk '{print $3, $6/($1*$1*$1)}' $file > CP_input
  awk '{print $3, $6/($1*$1*$1), $7/($1*$1*$1*$1*$1*$1), $16}' $file > CP_input
  ./CPderivate
  mv CP_derivate.out L$size\_P$pres\_CP.dat
  rm CP_input
  #######
#  if [ -f  L$size\_P$pres\_Volume.dat.SI ]; then
#     awk '{print $1, $3*1E32/(1.38*$1*$2)}' L$size\_P$pres\_Volume.dat.SI > L$size\_P$pres\_KT.dat.SI  ## KT en 10^-3 MPa^-1  
#  fi
  ########################
#  file=L$size\_P$pres\_cluster.dat
#  if [ -f $file ]; then
#   awk '{print $3, $4/(6*$1*$1*$1)}' $file > L$size\_P$pres\_largest_cluster.dat
#     sort -g L$size\_P$pres\_largest_cluster.dat > L$size\_P$pres\_largest_cluster.aux
#     mv L$size\_P$pres\_largest_cluster.aux L$size\_P$pres\_largest_cluster.dat  
#   awk '{print $3, $6/(6*$1*$1*$1)}' $file > L$size\_P$pres\_second_cluster.dat
#     sort -g L$size\_P$pres\_second_cluster.dat > L$size\_P$pres\_second_cluster.aux
#     mv L$size\_P$pres\_second_cluster.aux L$size\_P$pres\_second_cluster.dat  
#   awk '{print $3, $8}' $file > L$size\_P$pres\_percolation_probability.dat
#     sort -g L$size\_P$pres\_percolation_probability.dat > L$size\_P$pres\_percolation_probability.aux
#     mv L$size\_P$pres\_percolation_probability.aux L$size\_P$pres\_percolation_probability.dat  
#   awk '{print $3, $9/(6*$1*$1*$1)}' $file > L$size\_P$pres\_average_cluster_size.dat
#     sort -g L$size\_P$pres\_average_cluster_size.dat > L$size\_P$pres\_average_cluster_size.aux
#     mv L$size\_P$pres\_average_cluster_size.aux L$size\_P$pres\_average_cluster_size.dat
#   awk '{print $3, $11/(6*$1*$1*$1)}' $file > L$size\_P$pres\_average_num_clusters.dat
#     sort -g L$size\_P$pres\_average_num_clusters.dat > L$size\_P$pres\_average_num_clusters.aux
#     mv L$size\_P$pres\_average_num_clusters.aux L$size\_P$pres\_average_num_clusters.dat   
#  fi
done
done
######
#for s1 in $L
#do
#if [ -f KT_list_isobars ]; then
#  rm KT_list_isobars
#fi
#awk '{print $3, $2, ($1*$1*$1)/$8, ($1*$1*$1)*sqrt(3*$9/$16)/($8*$8) }' L$s1\_average_data.dat > KT_input
#for p1 in $P
#do
#  echo $p1 >> KT_list_isobars
#done
#./KTderivate
#rm KT_input KT_list_isobars
#for file in $(ls P*KTderivate.dat)
#do
#  sort -g $file > L$s1\_$file
#  rm $file
#done
#done
######

