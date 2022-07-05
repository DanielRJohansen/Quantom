#!/bin/bash

echo "#### LIMA-Dynamics installer is now beginning ####"


default_dir=~/Desktop/LIMA/

echo "Using $default_dir as install directory"
#rm -rf "$default_dir"/


mkdir -p "$default_dir/Applications"
mkdir -p "$default_dir/Simulation"



######### Now make Quantom
Q_dir="$default_dir"/Applications/Quantom
mkdir -p "$Q_dir"/include
mkdir "$Q_dir"/src
mkdir "$Q_dir"/build

cp ./Quantom/Quantom/*.*h "$Q_dir"/include/
cp ./Quantom/Quantom/*.cu "$Q_dir"/src/
cp ./Quantom/Quantom/*.cpp "$Q_dir"/src/

cp ./LIMA/Applications/Quantom/* "$Q_dir"/



# Handle GLFW for Quantom. This is a quite bad way...
GLFW_CMAKER=~/Downloads/glfw-3.3.6/CMakeLists.txt
if [ ! -f "$GLFW_CMAKER" ]; then
	echo "glfw-3.3.6 not present in Downloads folder. Terminating install."
	exit 1
fi
mkdir -p "$default_dir"/Dependencies/GLFW/
cp -r ~/Downloads/glfw-3.3.6/* "$default_dir"/Dependencies/GLFW/


#exit 0

######### Now make FFM
FFM_dir="$default_dir"/Applications/Forcefieldmaker
mkdir -p "$FFM_dir"/include
mkdir "$FFM_dir"/src
mkdir "$FFM_dir"/build

cp ./LIMA_ForcefieldMaker/LIMA_ForcefieldMaker/*.h "$FFM_dir"/include/
cp ./LIMA_ForcefieldMaker/LIMA_ForcefieldMaker/*.cpp "$FFM_dir"/src/
cp ./LIMA/Applications/Forcefieldmaker/* "$FFM_dir"/				# CMakeLists.txt and build.sh


## Now make SimPostprocessor
SPP_dir="$default_dir"/Applications/SimPreprocessor
mkdir -p "$SPP_dir"/include
mkdir "$SPP_dir"/src
mkdir "$SPP_dir"/build

cp ./LIMA_services/LIMA_services/*.h "$SPP_dir"/include/
cp ./LIMA_services/LIMA_services/*.cpp "$SPP_dir"/src/


## Make test sim
s_dir="$default_dir"/Simulation
mkdir -p "$s_dir/"Molecule
mkdir "$s_dir"/Forcefield
#cp ~/Downloads/QnD/* "$s_dir"/Forcefield/	# nb and b ff
#cp ~/Downloads/QnD/*/* "$s_dir"/Molecule/	# conf and topol
cp ./Demo/6lzm/* "$s_dir"/Molecule/		# conf and topol"
cp ./Demo/CHARMM_FF/* "$s_dir"/Forcefield/	# nb and b ff



## Compile all applications
cd "$FFM_dir"
chmod +x build.sh
./build.sh
#./ffmrun




cd "$Q_dir"
chmod +x build.sh
chmod +x mdrun.sh
./build.sh

printf "\n\n\n\n\n\n\n\n\n\n"
printf "All LIMA applications have been installed\n"



## Run DEMO

#read -p "Press y to start demo simulation    " confirm && [[ $confirm == [yY] ]] || exit 1

cd "$FFM_dir"
./ffmrun


cd "$Q_dir"
./mdrun.sh



