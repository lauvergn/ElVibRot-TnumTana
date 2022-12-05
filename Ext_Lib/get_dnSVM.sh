rm -r AD_dnSVM*
rm -f dnSVMLib #always remove the link

ReleaseVersion=Save_AD_dnSVM-2.0.0.zip


#latest release
#latest HEAD version
 version=https://github.com/lauvergn/AD_dnSVM/archive/refs/heads/main.zip


curl -LJ $version --output dnSVM.zip
test -e dnSVM.zip && echo dnSVM.zip file exist || cp $ReleaseVersion dnSVM.zip
unzip dnSVM.zip
rm -f dnSVM.zip

QMLDIR=`ls -d AD_dnSVM*`
#echo $QMLDIR

ln -s $QMLDIR dnSVMLib
