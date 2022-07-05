#!/bin/bash


if [ ! -d ./build ];
then
	echo "Build folder not found"
	 return;
fi

cd build/

cmake ../
make

mv ./ffmrun ../
