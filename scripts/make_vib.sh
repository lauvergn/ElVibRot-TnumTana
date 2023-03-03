#!/bin/bash

DIR_EVRT=$1
F90=$2

DIR_pot=$DIR_EVRT/Ext_Lib/Coord_KEO_PrimOp/sub_pot
DIR_T=$DIR_EVRT/Ext_Lib/Coord_KEO_PrimOp/Source_TnumTana_Coord/sub_operator_T

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

ls $DIR_pot

if (diff $DIR_pot/sub_system_save.f $DIR_pot/sub_system.f >/dev/null)
then 
  echo "the sub_system.f files are identical"
else 
  cp $DIR_pot/sub_system_save.f   $DIR_pot/sub_system.f
  echo "the sub_system.f files are different"
fi
#cp $DIR_pot/sub_system_save.f   $DIR_pot/sub_system.f 

if (diff $DIR_pot/sub_system_save.f90 $DIR_pot/sub_system.f90 >/dev/null)
then 
  echo "the sub_system.f90 files are identical"
else 
  cp $DIR_pot/sub_system_save.f90 $DIR_pot/sub_system.f90
  echo "the sub_system.f90 files are different"
fi
#cp $DIR_pot/sub_system_save.f90 $DIR_pot/sub_system.f90

if (diff $DIR_T/calc_f2_f1Q_save.f90 $DIR_T/calc_f2_f1Q.f90 >/dev/null)
then 
  echo "the calc_f2_f1Q.f90 files are identical"
else 
  cp $DIR_T/calc_f2_f1Q_save.f90 $DIR_T/calc_f2_f1Q.f90
  echo "the calc_f2_f1Q.f90 files are different"
fi

if (diff $DIR_T/Calc_Tab_dnQflex_save.f90 $DIR_T/Calc_Tab_dnQflex.f90 >/dev/null)
then 
  echo "the Calc_Tab_dnQflex.f90 files are identical"
else 
  cp $DIR_T/Calc_Tab_dnQflex_save.f90 $DIR_T/Calc_Tab_dnQflex.f90
  echo "the Calc_Tab_dnQflex files are different"
fi

if (diff $DIR_T/Sub_X_TO_Q_ana_save.f90 $DIR_T/Sub_X_TO_Q_ana.f90 >/dev/null)
then 
  echo "the Sub_X_TO_Q_ana.f90 files are identical"
else 
  cp $DIR_T/Sub_X_TO_Q_ana_save.f90 $DIR_T/Sub_X_TO_Q_ana.f90
  echo "the Sub_X_TO_Q_ana.f90 files are different"
fi


rm \$name_file " >> vib