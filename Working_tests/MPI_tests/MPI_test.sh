#!/bin/bash

here=$(pwd)

## check OS and available processor
if [ $(uname) == 'Darwin' ]; then
  num_cores=$(sysctl -n hw.physicalcpu)
  echo 'test code MPI on MAC OS, available processors:' $num_cores
elif [ $(uname) == 'Linux' ]; then
  num_cores=$(grep -c ^processor /proc/cpuinfo)
  echo 'test code MPI on Linux, available processors:' $num_cores
fi

if [ $num_cores -lt  2 ]; then
  echo 'Warning: not enough processors for testing MPI'
  exit
fi

echo 'test log in MPI_test.log'

## -------------------------------------------------------------------------------------
echo 'Davidson example:'
## -------------------------------------------------------------------------------------

## Davidson test: 6D
## -------------------------------------------------------------------------------------
## 6D: scheme 0-2
for ii in {0..2}
do
  echo "> 6D, MPI scheme $ii: result in $here/6D_Davidson_S$ii/result"
  cd 6D_Davidson_S$ii

  ./run_jobs > MPI_test.log
  file="result/res_HenonHeiles_6D_SGtype4_LB3_B2_LG3_2MPIcores"
  if [ -f "$file" ]
  then

    echo "$file found."
    grep 'lev0' "$file" | awk '{print $5}' > ./result/levels
    err=$(awk '{if(NR==FNR){bench1[NR]=$1} else {if(bench1[FNR]-$1>0.00000001 || bench1[FNR]-$1<-0.00000001){print '1'}}}'  benchmark ./result/levels)
    if [ -z $err ]; then echo '  PASS'; else echo "ERROR in 6D Davidson MPI scheme $ii test"; fi

  else
    echo "ERROR in 6D Davidson MPI scheme $ii test, $file not found."
  fi

  cd ..
done

## 6D: scheme 3
if [ $num_cores -gt  2 ]; then
  echo "> 6D, MPI quasi-scheme 3: result in $here/6D_Davidson_S3/result"
  cd 6D_Davidson_S3

  ./run_jobs >> MPI_test.log
  file="result/res_HenonHeiles_6D_SGtype4_LB3_B2_LG3_3MPIcores"
  if [ -f "$file" ]
  then 

    echo "$file found."
    grep 'lev0' "$file" | awk '{print $5}' > ./result/levels
    err=$(awk '{if(NR==FNR){bench1[NR]=$1} else {if(bench1[FNR]-$1>0.00000001 || bench1[FNR]-$1<-0.00000001){print '1'}}}'  benchmark ./result/levels)
    if [ -z $err ]; then echo "  PASS"; else echo "ERROR in 6D Davidson MPI scheme 3 test"; fi

  else
    echo "ERROR in 6D Davidson MPI scheme 3 test, $file not found."
  fi

  cd ..

else
  echo 'not enough processors for a simple test of MPI scheme 3'
fi
## -------------------------------------------------------------------------------------


## Davidson test: 21D
## -------------------------------------------------------------------------------------
## 21D: scheme 0-2
for ii in {0..2}
do
  echo "> 21D, MPI scheme $ii: result in $here/21D_Davidson_S$ii/result"
  cd ./21D_Davidson_S$ii

  ./run_jobs >> MPI_test.log
  file="result/res_HenonHeiles_21D_SGtype4_LB2_B2_LG2_2MPIcores"
  if [ -f "$file" ]
  then

    echo "$file found."
    grep 'lev0' "$file" | awk '{print $5}' > ./result/levels
    err=$(awk '{if(NR==FNR){bench1[NR]=$1} else {if(bench1[FNR]-$1>0.00000001 || bench1[FNR]-$1<-0.00000001){print '1'}}}'  benchmark ./result/levels)
    if [ -z $err ]; then echo "  PASS"; else echo "ERROR in 21D Davidson MPI scheme $ii test"; fi

  else
    echo "ERROR in 21D Davidson MPI scheme $ii test, $file not found."
  fi

  cd ..
done

## 21D: scheme 3
if [ $num_cores -gt  2 ]; then
  echo "> 21D, MPI quasi-scheme 3: result in $here/21D_Davidson_S3/result"
  cd 21D_Davidson_S3

  ./run_jobs >> MPI_test.log
  file="result/res_HenonHeiles_21D_SGtype4_LB2_B2_LG2_3MPIcores"
  if [ -f "$file" ]
  then

    echo "$file found."
    grep 'lev0' "$file" | awk '{print $5}' > ./result/levels
    err=$(awk '{if(NR==FNR){bench1[NR]=$1} else {if(bench1[FNR]-$1>0.00000001 || bench1[FNR]-$1<-0.00000001){print '1'}}}'  benchmark ./result/levels)
    if [ -z $err ]; then echo "  PASS"; else echo "ERROR in 21D Davidson MPI scheme 3 test"; fi

  else
    echo "ERROR in 21D Davidson MPI scheme 3 test, $file not found."
  fi

  cd ..

else
  echo 'not enough processors for a simple test of MPI scheme 3'
fi
## -------------------------------------------------------------------------------------


## -------------------------------------------------------------------------------------
echo 'Propagation example:'
## -------------------------------------------------------------------------------------

## Propagation test: 12D
## -------------------------------------------------------------------------------------
# 12D: scheme 0-2
for ii in {0..2}
do
  echo "> 12D, MPI scheme $ii: result in $here/12D_propagation_S$ii/result"
  cd 12D_propagation_S$ii

  ./run_jobs >> MPI_test.log
  file="result/file_auto"
  if [ -f "$file" ]
  then

    echo "$file found."
    grep 'AutoCor' "$file" | awk '{print $2, $3, $4, $5}' > ./result/auto_cor
    err=$(awk '{if(NR==FNR){bench2[NR]=$2; bench3[NR]=$3; bench4[NR]=$4} else {if(bench2[FNR]-$2>0.00000001 || bench2[FNR]-$2<-0.00000001 || bench3[FNR]-$3>0.00000001 || bench3[FNR]-$3<-0.00000001 || bench4[FNR]-$4>0.00000001 || bench4[FNR]-$4<-0.00000001){print '1'}}}'  benchmark ./result/auto_cor)
    if [ -z $err ]; then echo "  PASS"; else echo "ERROR in 12D propagation MPI scheme $ii test"; fi

  else
    echo "ERROR in 12D propagation MPI scheme $ii test, $file not found."
  fi

  cd ..
done

# 12D: scheme 3
if [ $num_cores -gt  2 ]; then
  echo "> 12D, MPI quasi-scheme 3: result in $here/12D_propagation_S3/result"
  cd 12D_propagation_S3

  ./run_jobs >> MPI_test.log
  file="result/file_auto"
  if [ -f "$file" ]
  then

    echo "$file found."
    grep 'AutoCor' "$file" | awk '{print $2, $3, $4, $5}' > ./result/auto_cor
    err=$(awk '{if(NR==FNR){bench2[NR]=$2; bench3[NR]=$3; bench4[NR]=$4} else {if(bench2[FNR]-$2>0.00000001 || bench2[FNR]-$2<-0.00000001 || bench3[FNR]-$3>0.00000001 || bench3[FNR]-$3<-0.00000001 || bench4[FNR]-$4>0.00000001 || bench4[FNR]-$4<-0.00000001){print '1'}}}'  benchmark ./result/auto_cor)
    if [ -z $err ]; then echo "  PASS"; else echo "ERROR in 12D propagation MPI scheme 3 test"; fi

  else
    echo "ERROR in 12D propagation MPI scheme 3 test, $file not found."
  fi

  cd ..

else
  echo 'not enough processors for a simple test of MPI scheme 3'
fi
## -------------------------------------------------------------------------------------


## Propagation test: 24D
## ----------------------------------------------------------------------------------
# 24D: scheme 0-2
for ii in {0..2}
do
  echo "> 24D, MPI scheme $ii: result in $here/24D_propagation_S$ii/result"
  cd 24D_propagation_S$ii

  ./run_jobs >> MPI_test.log
  file="result/file_auto"
  if [ -f "$file" ]
  then

    echo "$file found."
    grep 'AutoCor' "$file" | awk '{print $2, $3, $4, $5}' > ./result/auto_cor
    err=$(awk '{if(NR==FNR){bench2[NR]=$2; bench3[NR]=$3; bench4[NR]=$4} else {if(bench2[FNR]-$2>0.00000001 || bench2[FNR]-$2<-0.00000001 || bench3[FNR]-$3>0.00000001 || bench3[FNR]-$3<-0.00000001 || bench4[FNR]-$4>0.00000001 || bench4[FNR]-$4<-0.00000001){print '1'}}}'  benchmark ./result/auto_cor)
    if [ -z $err ]; then echo "  PASS"; else echo "ERROR in 24D propagation MPI scheme $ii test"; fi

  else
    echo "ERROR in 24D propagation MPI scheme $ii test, $file not found."
  fi

  cd ..
done

# 24D: scheme 3
if [ $num_cores -gt  2 ]; then
  echo "> 24D, MPI quasi-scheme 3: result in $here/24D_propagation_S3/result"
  cd 24D_propagation_S3

  ./run_jobs >> MPI_test.log
  file="result/file_auto"
  if [ -f "$file" ]
  then

  echo "$file found."
  grep 'AutoCor' "$file" | awk '{print $2, $3, $4, $5}' > ./result/auto_cor
  err=$(awk '{if(NR==FNR){bench2[NR]=$2; bench3[NR]=$3; bench4[NR]=$4} else {if(bench2[FNR]-$2>0.00000001 || bench2[FNR]-$2<-0.00000001 || bench3[FNR]-$3>0.00000001 || bench3[FNR]-$3<-0.00000001 || bench4[FNR]-$4>0.00000001 || bench4[FNR]-$4<-0.00000001){print '1'}}}'  benchmark ./result/auto_cor)
  if [ -z $err ]; then echo '  PASS'; else echo "ERROR in 24D propagation MPI scheme 3 test"; fi

  else
    echo "ERROR in 24D propagation MPI scheme 3 test, $file not found."
  fi

  cd ..

else
  echo "not enough processors for a simple test of MPI scheme 3"
fi


## ----------------------------------------------------------------------------------
## scheme 3
echo "To test MPI scheme 3, set MPI_fake_nodes>0 with more processores or submit jobs to more than one node"











