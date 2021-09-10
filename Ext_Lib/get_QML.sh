rm -r QuantumModelLib*
rm -f QuantumModelLib #always remove the link

#latest release
 version=https://github.com/lauvergn/QuantumModelLib/archive/refs/tags/v7.3.zip
#latest HEAD version
#version=https://github.com/lauvergn/QuantumModelLib/archive/refs/heads/OOP_branch.zip


curl -LJ $version --output OOP_branch.zip
unzip OOP_branch.zip
rm -f OOP_branch.zip

QMLDIR=`ls -d QuantumModelLib*`
#echo $QMLDIR

ln -s $QMLDIR QuantumModelLib

