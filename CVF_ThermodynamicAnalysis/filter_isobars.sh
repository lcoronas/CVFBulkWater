# This scripts filters data in L32_average_data.dat to generate L32_P*.dat file for every isobar 

L="32"
P="-0.8 -0.75 -0.7 -0.6 -0.5 -0.4 -0.38 -0.35 -0.3 -0.25 -0.22 -0.2 -0.15 -0.13 -0.1 -0.05 0.0 0.05 0.1 0.13 0.15 0.2 0.24 0.25 0.3 0.343 0.35 0.38 0.4 0.45 0.47 0.5 0.515 0.53 0.55 0.6 0.63 0.66 0.7 0.73 0.77 0.8 0.85 0.88 0.9 0.95 1"
for size in $L
do
file=L$size\_average_data.dat
sed "s|\t|  |g" $file > average_data_aux.dat
cp average_data_aux.dat $file
rm average_data_aux.dat
###########
for pres in $P
do
    echo "Isobar L"$size" P"$pres"..."
    grep "^$size  $pres " $file >> P$pres
    sort -g P$pres > L$size\_P$pres.dat
    rm P$pres
done
done
