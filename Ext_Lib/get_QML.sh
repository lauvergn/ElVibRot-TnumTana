#!/bin/bash

EXTLIB_TYPE=$1

echo "In get_QML.sh"


SAVE_version=Save_QuantumModelLib-20.1-dev
LOC_version=QuantumModelLib

rm -r QuantumModelLib*
rm -f QuantumModelLib #always remove the link



#latest release
#latest HEAD version (dev version)
 version=https://github.com/lauvergn/QuantumModelLib/archive/refs/tags/v20.1-dev.zip


test -z $EXTLIB_TYPE       &&    curl -LJ $version --output $LOC_version.zip
test $EXTLIB_TYPE != 'loc' &&    curl -LJ $version --output $LOC_version.zip

test -e $LOC_version.zip && echo $LOC_version.zip file exist || cp $SAVE_version.zip $LOC_version.zip

unzip $LOC_version.zip
rm -f $LOC_version.zip



QMLDIR=`ls -d QuantumModelLib*`
#echo $QMLDIR

ln -s $QMLDIR QuantumModelLib

echo "End get_QML.sh"
