#!/bin/sh


if [ ! -d ./build ];
then
	echo "Build directory not found"
	return
fi

cd build
rm -rf ./*



cmake ../
make

mv mdrun ../

#./lima
#nohup ./lima &


