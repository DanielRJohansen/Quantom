#!/bin/bash


ff=../../Simulation/Forcefield/Forcefield.txt
ffs=../../Simulation/Forcefield/ForcefieldSummary.txt

conf=../../Simulation/Molecule/conf.gro
topol=../../Simulation/Molecule/topol.top


if [[ ! -f "$conf" || ! -f "$topol" ]];
then
	echo "Conf or Topol file is missing"
	exit 1
fi


if [[ ! -f "$ff" || !  -f "$ffs" ]];
then
	echo "Making forcefield files now"
	../Forcefieldmaker/ffm
fi

./mdrun
