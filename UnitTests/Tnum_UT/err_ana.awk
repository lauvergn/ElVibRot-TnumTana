BEGIN {err=0} 
{
  for (i=2;i<=NF;i++) {err=err+$i*$i}
} 
END {
  err=sqrt(err)
  if (err < 1.e-8) {print "NO ERROR in " $1} 
  else {print $0 ; print "ERROR(s) in " $1}
}
