#!/bin/bash

file=$1
epsi=$2
col=$3
info=$4

 LANG=C awk -v epsi=$epsi -v col=$col -v info=$info '
   BEGIN {x=0}
         {v=$col ; x=v-x}
   END {
     e=sqrt(x*x)
     print "For " info ", Val: " v " , Diff: " e
     if (e > epsi)
       {print "    ERROR for " info ",  largest diff > " epsi}
     else
       {print "    No PROBLEM for " info ", largest diff <= " epsi}
  }' $file
