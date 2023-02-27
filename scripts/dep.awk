#MODULE mod_auto_basis sub_Basis/sub_Auto_Basis.f90
BEGIN {
  mod=0
} 
{
  if (tolower($1) == "use" && tolower($2) == "mod_auto_basis") {
     mod = 1
  }
}
END {
  print "file: " FILENAME
  n=split(FILENAME,tab,"/")
  if (n > 0) {
    print n " " tab[n]
    l=length(tab[n])
    objf=substr(tab[n],1,l-4)
  }
  if (mod == 1)  print "$(OBJ_DIR)/" objf ".o : $(OBJ_DIR)/sub_Auto_Basis.o"
}
