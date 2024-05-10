# PI = 3.1415927
# PI/2 = 1.5707693
# PI/4 = 0.7853982
############
L="8"
###
#T2="0.0166 0.0167 0.0168 0.0169 0.017 0.0171 0.0172 0.0173"
#P2="0.8"
###
#T2="0.065 0.0651 0.0652 0.0653 0.0654 0.0655 0.0656 0.0657 0.0658 0.0659 0.066 0.0661 0.0662 0.0663 0.0664 0.0665 0.0666 0.0667 0.0668 0.0669"
#P2="0.7"
###
#T2="0.0924 0.0925 0.0926 0.0927 0.0928 0.0929 0.093 0.0931 0.0932 0.0933 0.0934 0.0935"
#P2="0.6"
###
T2="0.1058 0.1059 0.106 0.1061 0.1062 0.1063 0.1064 0.1065 0.1066 0.1067 0.1068"
P2="0.5"
###
#T2="0.1152 0.1153 0.1154 0.1155 0.1156 0.1157 0.1158 0.1159"
#P2="0.4"
###
#T2="0.127 0.1271 0.1272 0.1273 0.1274 0.1275 0.1276 0.1277 0.1278 0.1279"
#P2="0.2"
###
#T2="0.134 0.1341 0.1342 0.1343 0.1344 0.1345 0.1346 0.1347 0.1348 0.1349 0.135"
#P2="0.0"
###
#T2="0.1376 0.1377 0.1378 0.1379 0.138 0.1381 0.1382 0.1383 0.1384 0.1385 0.1386 0.1387 0.1388"
#P2="-0.2"
###
#T2="0.1404 0.1405 0.1406 0.1407 0.1408 0.1409 0.141 0.1411 0.1412 0.1413 0.1414 0.1415 0.1416 0.1417 0.1418 0.1419 0.142 0.1421 0.1422"
#P2="-0.7"
############
#L="12"
###
#T2="0.0166 0.0167 0.0168 0.0169 0.017"
#P2="0.8"
###
#T2="0.0652 0.0653 0.0654 0.0655 0.0656 0.0657 0.0658"
#P2="0.7"
###
#T2="0.0818 0.0819 0.082 0.0821 0.0822"
#P2="0.65"
###
#T2="0.0924 0.0925 0.0926 0.0927 0.0928"
#P2="0.6"
###
#T4="0.0999 0.1 0.1001 0.1002 0.1003"
#P4="0.55"
###
#T3="0.106 0.1061 0.1062"
#P3="0.5"
###
#T2="0.1109 0.111 0.1111 0.1112"
#P2="0.45"
###
#T4="0.113 0.1131 0.1132 0.1133 0.1134 0.1154"
#T4="0.1154"
#P4="0.425"
###
#T2="0.1151 0.1152 0.1153 0.1154"
#P2="0.4"
###
#T3="0.115 0.1169 0.117 0.1171 0.1172 0.1173"
#T3="0.115"
#P3="0.375"
###
#T2="0.1339 0.134 0.1341 0.1342"
#P2="0.0"
###
#T2="0.1406 0.1407 0.1408 0.1409 0.141 0.1411"
#P2="-0.7"
############
#L="16"
###
#T2="0.0166 0.0167 0.0168"
#P2="0.8"
###
#T2="0.0418 0.0419 0.042"
#P2="0.75"
###
#T2="0.0564 0.0565 0.0566"
#P2="0.72"
###
#T2="0.064 0.0651 0.0652 0.0653 0.0654 0.0655 0.068"
#T2="0.0651 0.0652 0.0653 0.0654 0.0655"
#P2="0.7"
###
#T2="0.0728 0.0729 0.0730"
#P2="0.68"
###
############
N_SIMULATION="11"
N_MOLECULES=$(($L*$L*$L))
#Temp: 8 decimales maximo. DeltaT: 10 dec max
#MIN_T="0.1040"
#MAX_T="0.1080"
###### P=0.8
#MIN_T="0.0150"
#MAX_T="0.0190"
###### P=0.75
#MIN_T="0.04185"
#MAX_T="0.04191"
###### P=0.72
#MIN_T="0.056456"
#MAX_T="0.056600"
###### P=0.68
#MIN_T="0.07288"
#MAX_T="0.07297"
###### P=0.7
#MIN_T="0.065288"
#MAX_T="0.066000"
###### P=0.65
#MIN_T="0.08201"
#MAX_T="0.08210"
###### P=0.6
#MIN_T="0.092575"
#MAX_T="0.09300"
###### P=0.5
#MIN_T="0.106099"
#MAX_T="0.106400"
###### P=0.5
MIN_T="0.106380"
MAX_T="0.106450"
###### P=0.4
#MIN_T="0.115525"
#MAX_T="0.116000"
#MIN_T="0.115207"
#MAX_T="0.115300"
###### P=0.2
#MIN_T="0.127485"
#MAX_T="0.127700"
###### P=0.0
#MIN_T="0.13445"
#MAX_T="0.13490"
###### P=-0.2
#MIN_T="0.13824"
#MAX_T="0.13880"
###### P=-0.7
#MIN_T="0.1407"
#MAX_T="0.1418"
DELTA_T="0.2"
#Pres: 5 decimales maximo. DeltaP: 8 dec max
MIN_P="0.500"
MAX_P="0.501"
DELTA_P="0.1"
#DELTA_P="0.1"
#Mixing Parameter: 5 decimales maximo. DeltaS: 7 dec max
MIN_S="1.588"
MAX_S="1.589"
######  mixing parameter
#MIN_S="1.6280"    #L=8 P=0.2 T=0.127485 s=1.628
#MAX_S="1.6281"
#MIN_S="1.5845"   #L=12 P=0.6 T=0.092575 s=1.5845
#MAX_S="1.6000"
#MIN_S="1.57775"  #L=16 P=0.7 T=0.065288 s=1.57775
#MAX_S="1.60000"
DELTA_S="0.1"
###### alpha
#MIN_S="5.9000"
#MAX_S="6.1000" 
#DELTA_S="0.01"
#### BIN SIZE FACTOR   n_bin2D = (int) ( binFactor * 10 * sqrt ( ndatos ) / 2 )
####                   n_bin_mixed = sqrt( n_bin2D )  
#BIN_FCT="0.1"
#BIN_FCT="1.0"
BIN_FCT="2.0"
#################
################# CALCULATE BINDER? 0 NO; 1 YES
#binder=1
echo "N_SIMULACIONES" $N_SIMULATION > input
echo "STEP 300" >> input
echo "L" $L >> input
echo "MIN_T" $MIN_T >> input
echo "MAX_T" $MAX_T >> input
echo "DELTA_T" $DELTA_T >> input
echo "MIN_P" $MIN_P >> input
echo "MAX_P" $MAX_P >> input
echo "DELTA_P" $DELTA_P >> input
echo "MIN_S" $MIN_S >> input
echo "MAX_S" $MAX_S >> input
echo "DELTA_S" $DELTA_S >> input
echo "BIN_FCT" $BIN_FCT >> input
let counter=0
for temp in $T2
do
  temp1=$temp
  x=$( echo -n $temp1 | wc -c )
  while [ !  $x == 9  ]
  do
    temp1=$( echo $temp1\0 )
    x=$( echo -n $temp1 | wc -c )
  done 
  for pres in $P2
  do
    if [ -f L$L\_T$temp\_P$pres/data.out ]; then
      pres1=$pres
      y=$( echo -n $pres1 | wc -c )
      while [ !  $y == 6  ]
      do
        pres1=$( echo $pres1\0 )
        y=$( echo -n $pres1 | wc -c )
      done
      cp L$L\_T$temp\_P$pres/data.out .
      mv data.out data_L$L\_T$temp1\_P$pres1.aux
      awk -v N=$N_MOLECULES '{print $3/N, N/$4}' data_L$L\_T$temp1\_P$pres1.aux > data_L$L\_T$temp1\_P$pres1
      rm data_L$L\_T$temp1\_P$pres1.aux
      NUMOFLINES=$(wc -l < data_L$L\_T$temp1\_P$pres1)
      echo $temp1 $pres1 $NUMOFLINES >> input
      let counter=$counter+1
    fi
  done
done
for temp in $T3
do
  temp1=$temp
  x=$( echo -n $temp1 | wc -c )
  while [ !  $x == 9  ]
  do
    temp1=$( echo $temp1\0 )
    x=$( echo -n $temp1 | wc -c )
  done 
  for pres in $P3
  do
    if [ -f L$L\_T$temp\_P$pres/data.out ]; then
      pres1=$pres
      y=$( echo -n $pres1 | wc -c )
      while [ !  $y == 6  ]
      do
        pres1=$( echo $pres1\0 )
        y=$( echo -n $pres1 | wc -c )
      done
      cp L$L\_T$temp\_P$pres/data.out .
      mv data.out data_L$L\_T$temp1\_P$pres1.aux
      awk -v N=$N_MOLECULES '{print $3/N, N/$4}' data_L$L\_T$temp1\_P$pres1.aux > data_L$L\_T$temp1\_P$pres1
      rm data_L$L\_T$temp1\_P$pres1.aux
      NUMOFLINES=$(wc -l < data_L$L\_T$temp1\_P$pres1)
      echo $temp1 $pres1 $NUMOFLINES >> input
      let counter=$counter+1
    fi
  done
done
for temp in $T4
do
  temp1=$temp
  x=$( echo -n $temp1 | wc -c )
  while [ !  $x == 9  ]
  do
    temp1=$( echo $temp1\0 )
    x=$( echo -n $temp1 | wc -c )
  done 
  for pres in $P4
  do
    if [ -f L$L\_T$temp\_P$pres/data.out ]; then
      pres1=$pres
      y=$( echo -n $pres1 | wc -c )
      while [ !  $y == 6  ]
      do
        pres1=$( echo $pres1\0 )
        y=$( echo -n $pres1 | wc -c )
      done
      cp L$L\_T$temp\_P$pres/data.out .
      mv data.out data_L$L\_T$temp1\_P$pres1.aux
      awk -v N=$N_MOLECULES '{print $3/N, N/$4}' data_L$L\_T$temp1\_P$pres1.aux > data_L$L\_T$temp1\_P$pres1
      rm data_L$L\_T$temp1\_P$pres1.aux
      NUMOFLINES=$(wc -l < data_L$L\_T$temp1\_P$pres1)
      echo $temp1 $pres1 $NUMOFLINES >> input
      let counter=$counter+1
    fi
  done
done
#cat input
#./Print_histogram_reweighting_mixed_param
echo "Contador de simulaciones:" $counter
#echo "Starting hour:"
#date
./histogram_reweighting_mixed_paramSI_norm_x
#echo "Finishing hour:"
#date
#echo "ejecuto segundo programa"
#./histogram_reweighting_rotation
