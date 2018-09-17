#!/bin/bash
#
mkdir temp
cd temp
rm *
#~/bin/$ARCH/f77split ../quadpack.f
#
for FILE in `ls -1 *.f`;
do
  gfortran -c $FILE
  if [ $? -ne 0 ]; then
    echo "Errors compiling " $FILE
    exit
  fi
done
rm *.f
#
ar qc libquadpack.a *.o
rm *.o
#
mv libquadpack.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/$ARCH/libquadpack.a."
