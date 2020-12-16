#!/bin/bash

here=$(pwd)

## check OS and available processor
if [ $(uname) == 'Darwin' ]; then
  num_cores=$(sysctl -n hw.physicalcpu)
  echo 'test code openMP on MAC OS, available processors:' $num_cores
elif [ $(uname) == 'Linux' ]; then
  num_cores=$(grep -c ^processor /proc/cpuinfo)
  echo 'test code openMP on Linux, available processors:' $num_cores
fi

echo 'test log in openMP_test.log'

## ----------------------------------------------------------------------------------
echo 'Davidson test:'
## ----------------------------------------------------------------------------------

## Davidson test: 6D
## ----------------------------------------------------------------------------------
echo '> 6D, result in' $here'/6D_Davidson_openMP/result'
cd 6D_Davidson_openMP

./run_jobs > openMP_test.log
file="result/res_HenonHeiles_6D_SGtype4_LB3_B2_LG3_2openMpcores"
if [ -f "$file" ]
then

  echo "$file found."
  grep 'lev0' "$file" | awk '{print $5}' > ./result/levels
  err=$(awk '{if(NR==FNR){bench1[NR]=$1} else {if(bench1[FNR]-$1>0.00000001 || bench1[FNR]-$1<-0.00000001){print '1'}}}'  benchmark ./result/levels)
  if [ -z $err ]; then echo '  PASS'; else echo 'ERROR in 6D Davidson openMP test'; fi

else
	echo "ERROR in 6D Davidson openMP test, $file not found."
fi

cd ..



## Davidson test: 21D
## ----------------------------------------------------------------------------------
echo '> 21D, result in' $here'/21D_Davidson_openMP/result'
cd ./21D_Davidson_openMP

./run_jobs >> MPI_test.log
file="result/res_HenonHeiles_21D_SGtype4_LB2_B2_LG2_2openMpcores"
if [ -f "$file" ]
then

  echo "$file found."
  grep 'lev0' "$file" | awk '{print $5}' > ./result/levels
  err=$(awk '{if(NR==FNR){bench1[NR]=$1} else {if(bench1[FNR]-$1>0.00000001 || bench1[FNR]-$1<-0.00000001){print '1'}}}'  benchmark ./result/levels)
  if [ -z $err ]; then echo '  PASS'; else echo 'ERROR in 21D Davidson openMP test'; fi

else
	echo "ERROR in 21D Davidson openMP test, $file not found."
fi

cd ..



## ----------------------------------------------------------------------------------
echo 'Propagation test:'
## ----------------------------------------------------------------------------------

## Propagation test: 12D
## ----------------------------------------------------------------------------------
# 12D: scheme 1
echo '> 12D, result in' $here'/12D_Davidson_openMP/result'
cd 12D_propagation_openMP

./run_jobs >> MPI_test.log
file="result/file_auto"
if [ -f "$file" ]
then

  echo "$file found."
  grep 'AutoCor' "$file" | awk '{print $2, $3, $4, $5}' > ./result/auto_cor
  err=$(awk '{if(NR==FNR){bench2[NR]=$2; bench3[NR]=$3; bench4[NR]=$4} else {if(bench2[FNR]-$2>0.00000001 || bench2[FNR]-$2<-0.00000001 || bench3[FNR]-$3>0.00000001 || bench3[FNR]-$3<-0.00000001 || bench4[FNR]-$4>0.00000001 || bench4[FNR]-$4<-0.00000001){print '1'}}}'  benchmark ./result/auto_cor)
  if [ -z $err ]; then echo '  PASS'; else echo 'ERROR in 12D propagation openMP test'; fi

else
	echo "ERROR in 12D propagation openMP test, $file not found."
fi

cd ..

## Propagation test: 24D
## ----------------------------------------------------------------------------------
# 24D: scheme 1
echo '> 24D, result in' $here'/24D_propagation_openMP/result'
cd 24D_propagation_openMP

./run_jobs >> MPI_test.log
file="result/file_auto"
if [ -f "$file" ]
then

  echo "$file found."
  grep 'AutoCor' "$file" | awk '{print $2, $3, $4, $5}' > ./result/auto_cor
  err=$(awk '{if(NR==FNR){bench2[NR]=$2; bench3[NR]=$3; bench4[NR]=$4} else {if(bench2[FNR]-$2>0.00000001 || bench2[FNR]-$2<-0.00000001 || bench3[FNR]-$3>0.00000001 || bench3[FNR]-$3<-0.00000001 || bench4[FNR]-$4>0.00000001 || bench4[FNR]-$4<-0.00000001){print '1'}}}'  benchmark ./result/auto_cor)
  if [ -z $err ]; then echo '  PASS'; else echo 'ERROR in 24D propagation openMP test'; fi

else
	echo "ERROR in 24D propagation openMP test, $file not found."
fi

cd ..









