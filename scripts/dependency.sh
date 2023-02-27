#!/bin/bash

name_dep=1
SRCFile=2

list=`ls Source_ElVibRot/*.f90 Source_ElVibRot/*/*.f90 Source_ElVibRot/*/*/*.f90`

echo "#===============================================" > $name_dep
echo "#===============================================" > $SRCFile
echo "SRCFILE= \\" >> $SRCFile

for ff90 in $list
do
   awk -f scripts/mod2file.awk $ff90 >> $name_dep
   echo $ff90 | awk '{name=$1
   n=split(name,tab,"/")
  if (n > 0) {
    l=length(tab[n])
    print tab[n] " \\"
  }
  }' >> $SRCFile
done
echo "#===============================================" >> $name_dep
for ff90 in $list
do
   awk -f scripts/dep2.awk $ff90 >> $name_dep
done