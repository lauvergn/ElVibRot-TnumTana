#!/bin/bash

DIR_EVRT=$1
F90=$2

echo "#!/bin/bash
" > vib

if test $F90 = "pgf90"
then 
  echo "export OMP_STACKSIZE=50M" >> vib
else
  echo "#echo not pgf90" >> vib
fi

echo "
name_file='namelist'
i=0

while [ -e \$name_file ]
do
 i=\$(( \$i + 1 ))
 name_file='namelist'\$i
done
#echo \$name_file

cat > \$name_file

nice $DIR_EVRT/vib.exe --input \$name_file

rm \$name_file " >> vib
