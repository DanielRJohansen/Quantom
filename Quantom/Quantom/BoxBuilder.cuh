#pragma once
#include "Bodies.cuh"
#include "Simulation.cuh"


class BoxBuilder
{
public:
	BoxBuilder(Simulation* simulation) {

	}

private:
	int solvateBox();					// Returns # of solvate compounds placed
	Compound_H2O createCompound(Float3 com);
	bool spaceAvailable(Float3 com, float radius);
	void placeCompound(Compound_H2O);	// Simply assigns compounds to block - implement later
	void prepareCudaScheduler();

	Simulation* simToDevice();



	

	// ---------------------------------------------------- Variables ---------------------------------------------------- //
	Box box;




	// ---------------------------------------------------- Helper functions ---------------------------------------------------- //

};

