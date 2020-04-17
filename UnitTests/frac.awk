BEGIN{nb_err=0} 
{ if (NF > 3) {nb_err += ($(NF-2) != $NF)} } 
END {
if (nb_err == 0) 
   {print "NO ERROR in Frac UT"} 
else
   {print nb_err " ERRORS in Frac UT."; print "UnitTests/Check res_UT_Frac file."}}
