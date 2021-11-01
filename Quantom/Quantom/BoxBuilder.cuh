#pragma once
#include "Bodies.cuh"
#include "Simulation.cuh"


class BoxBuilder
{
public:
	BoxBuilder() {};
	void build(Simulation* simulation);

private:
	int solvateBox(Simulation* simulation);					// Returns # of solvate compounds placed
	Compound_H2O createCompound(Float3 com, int compound_index, CompoundState* statebuffer_node, CompoundNeighborList* neighborinfo_node);
	bool spaceAvailable(Float3 com, float radius);
	void compoundLinker(Simulation* simulation);									// Temp function


	

	// ---------------------------------------------------- Variables ---------------------------------------------------- //
	Box box;	// Local host version
	float box_len = BOX_LEN;
	float box_base = 0;
	const int max_compounds = 1'000'000;								// DO we want something const here? Smaller val?




	// We cannot use the pointers in the box, as they must be on device from start, 
	// since the compounds must know the adresses as they are created.
	CompoundState* compoundstates_host;		
	CompoundNeighborList* compoundneighborlists_host;
	//----------------------------------------------------------------------//
	
	



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

