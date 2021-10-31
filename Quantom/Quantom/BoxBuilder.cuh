#pragma once
#include "Bodies.cuh"
#include "Simulation.cuh"


class BoxBuilder
{
public:
	BoxBuilder(Simulation* simulation);

private:
	int solvateBox();					// Returns # of solvate compounds placed
	Compound_H2O createCompound(Float3 com, int compound_index, CompoundState* statebuffer_node, CompoundNeighborInfo* neighborinfo_node);
	bool spaceAvailable(Float3 com, float radius);
	void placeCompound(Compound_H2O compound);	// Simply assigns compounds to block - implement later
	void prepareCudaScheduler();

	Simulation* simToDevice();



	

	// ---------------------------------------------------- Variables ---------------------------------------------------- //
	Box box;	// Local host version
	float box_len = BOX_LEN;
	float box_base = 0;
	Simulation* simulation;	//Pointer to actual simulation



	// ---------------------------------------------------- Helper functions ---------------------------------------------------- //
	Float3 get3Random(int resolution) {
		return Float3(
			rand() % resolution / (float)resolution - 0.5,
			rand() % resolution / (float)resolution - 0.5,
			rand() % resolution / (float)resolution - 0.5
		);
	}
	float random() {
		return rand() % 10000 / 10000.f - 0.5f;
	}
};

