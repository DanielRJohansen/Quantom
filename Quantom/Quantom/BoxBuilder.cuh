#pragma once
#include "Bodies.cuh"
#include "Simulation.cuh"


class BoxBuilder
{
public:
	BoxBuilder() {};
	void build(Simulation* simulation);


private:

	void placeMainMolecule(Simulation* simulation);
	int solvateBox(Simulation* simulation);					// Returns # of solvate compounds placed
	Compound createCompound(Float3 com, int compound_index,
		CompoundState* statebuffer_node, CompoundNeighborList* neighborinfo_node, float dt);
	Solvent createSolvent(Float3 com, float dt);
	bool spaceAvailable(Box* box, Float3 com, float radius);
	void compoundLinker(Simulation* simulation);									// Temp function
	void solvateLinker(Simulation* simulation);
	void solvateCompoundCrosslinker(Simulation* simulation);
	

	// ---------------------------------------------------- Variables ---------------------------------------------------- //
	Box box;	// Local host version
	float box_len = BOX_LEN;
	float box_base = 0;





	float m = 18.01528;					// g/mol
	float M = m * 0.001;				// kg/mol
	//float k_B = 8.617333262145 * 10e-5;	// Boltzmann constant
	double k_B = 1.380 * 1e-23;
	float T = 293;	// Kelvin
	float R = 8.3144;					// J/(Kelvin*mol)
	float mean_velocity = m / (2 * k_B * T);
	float v_rms = sqrt(3 * R * T / M);










	// We cannot use the pointers in the box, as they must be on device from start, 
	// since the compounds must know the adresses as they are created.
	//CompoundState* compoundstates_host;		
	//CompoundNeighborList* compoundneighborlists_host;
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
		return rand() % 10000 / 10000.f * 2 - 1.f;
	}
};

