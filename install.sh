#!/bin/bash

echo "#### LIMA-Dynamics installer is now beginning ####"


default_dir=~/Desktop/LIMA/

echo "Using $default_dir as install directory"

mkdir -p "$default_dir/Applications"
mkdir -p "$default_dir/Simulation"



## Now make Quantom
Q_dir="$default_dir"/Applications/Quantom
mkdir -p "$Q_dir"/includes
mkdir "$Q_dir"/src
mkdir "Q_dir"/build

cp ./Quantom/Quantom/*.*h "$Q_dir"/includes/
cp ./Quantom/Quantom/*.cu "$Q_dir"/src/
cp ./Quantom/Quantom/*.cpp "$Q_dir"/src/

# Handle GLFW for Quantom. This is a quite bad way...
mkdir -p "$default_dir"/Dependencies/GLFW/
cp -r ~/Downloads/glfw-3.3.6/* "$default_dir"/Dependencies/GLFW/

## Now make FFM
FFM_dir="$default_dir"/Applications/Forcefieldmaker
mkdir -p "$FFM_dir"/includes
mkdir "$FFM_dir"/src
mkdir "$FFM_dir"/build

cp ./LIMA_ForcefieldMaker/LIMA_ForcefieldMaker/*.h "$FFM_dir"/includes/
cp ./LIMA_ForcefieldMaker/LIMA_ForcefieldMaker/*.cpp "$FFM_dir"/src/


## Now make SimPostprocessor
SPP_dir="$default_dir"/Applications/SimPreprocessor
mkdir -p "$SPP_dir"/includes
mkdir "$SPP_dir"/src
mkdir "$SPP_dir"/build

cp ./LIMA_services/LIMA_services/*.h "$SPP_dir"/includes/
cp ./LIMA_services/LIMA_services/*.cpp "$SPP_dir"/src/


