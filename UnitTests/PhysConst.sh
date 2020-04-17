res_UT=$1

nb_err=`grep -c "ERROR" $res_UT`
if [ $nb_err -eq 0 ]
then
  echo "NO ERROR in PhysConst UT"
else
  echo $nb_err " ERROR(S) in PhysConst UT"
  echo  " check the Examples/exa_PhysicalConstants/res* output files."
fi
