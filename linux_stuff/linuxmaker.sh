#!/bin/#!/bin/sh


mkdir -p ~/Desktop/LIMA/Dependencies/GLFW
mkdir ~/Desktop/LIMA/src
mkdir ~/Desktop/LIMA/include
mkdir ~/Desktop/LIMA/build

cp ./Quantom/Quantom/*.*h ~/Desktop/LIMA/include/
cp ./Quantom/Quantom/*.cpp ~/Desktop/LIMA/src/
cp ./Quantom/Quantom/*.cu ~/Desktop/LIMA/src/

cp -r ~/Downloads/glfw-3.3.6/* ~/desktop/LIMA/Dependencies/GLFW/
