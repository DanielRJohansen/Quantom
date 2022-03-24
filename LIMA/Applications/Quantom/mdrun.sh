#!/bin/bash


ff=../../Simulation/Forcefield.txt
ffs=../../Simulation/ForcefieldSummary.txt

conf=../../Simulation/conf.gro
topol=../../Simulation/topol.top

if [! [-f "$ff"] && [-f"$ffs"]];
then
	if [! [-f "$conf"] && [-f "$topol"]];
	then
		echo "Failed to make forcefield, missing conf.gro or topol.top file."
		return
	fi

	echo "Making forcefield files now"
	../Forcefieldmaker/ffm
fi

./mdrun
