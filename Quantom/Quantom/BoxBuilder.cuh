#pragma once
#include "Bodies.cuh"
#include "Simulation.cuh"
#include <vector>

class BoxBuilder
{
public:
	BoxBuilder() {
		srand(290128309);
	};
	void buildBox(Simulation* simulation);
	void addSingleMolecule(Simulation* simulation, Molecule* molecule);		// Can only use a single "add" function per Simulation for now!!!!!!!!!!!!!
	void addScatteredMolecules(Simulation* simulation, Compound* molecule, int n_copies);
	void addDoubleMembrane(Simulation* simulation, Compound* molecule);
	void finishBox(Simulation* simulation);
	int solvateBox(Simulation* simulation);					// Returns # of solvate compounds placed
	int solvateBox(Simulation* simulation, vector<Float3> *solvate_positions);	// Returns # of solvate compounds placed



private:

	void placeSingleMolecule(Simulation* simulation, Compound* molecule);

	void integrateCompound(Float3 com, int compound_index,
		CompoundState* statebuffer_node, double dt, Simulation* simulation);
	void integrateCompound(Compound* compound, Simulation* simulation);
	Solvent createSolvent(Float3 com, double dt);

	bool spaceAvailable(Box* box, Float3 com, double radius);
	void compoundLinker(Simulation* simulation);									// Temp function
	void solvateLinker(Simulation* simulation);
	void solvateCompoundCrosslinker(Simulation* simulation);
	


	// -------------- Functions for compound manipulation BEFORE integration -------------- //
	void placeMultipleCompoundsRandomly(Simulation* simulation, Compound* template_compound, int n_copies);
	Compound* randomizeCompound(Compound* template_compound);
	void moveCompound(Compound* compound, Float3 vector);

	//Float3 calcCompoundCom(Compound* compound);
	void rotateCompound(Compound* compound, Float3 xyz_rot);
	BoundingBox calcCompoundBoundingBox(Compound* compound);
	bool spaceAvailable(Box* box, Compound* compound);
	bool spaceAvailable(Box* box, Float3 particle_center);	// Ignore radius for this, as it just check against bounding boxes. 
	bool verifyPairwiseParticleMindist(Compound* a, Compound* b);
	//What about other solvents then? Not problem now while solvents are placed on a grid, but what about later?

	// ------------------------------------------------------------------------------------ //




	// ---------------------------------------------------- Variables ---------------------------------------------------- //
	Box box;	// Local host version
	double box_len = BOX_LEN;
	double box_base = 0;





	double M = SOLVENT_MASS;				// kg/mol
	//double k_B = 8.617333262145 * 10e-5;	// Boltzmann constant
	double k_B = 1.380 * 1e-23;
	double T = 313;	// Kelvin
	double R = 8.3144;					// J/(Kelvin*mol)
	//double mean_velocity = M / (2 * k_B * T);				// This is bullshit
	double v_rms = sqrt(3 * R * T / M);


	double MIN_NONBONDED_DIST = 0.5;



	Float3 most_recent_offset_applied = Float3(0.f);	// If molecule is offset, each solvent from .gro file must be aswell
														// Not a clean solution, i know.. :(



	// We cannot use the pointers in the box, as they must be on device from start, 
	// since the compounds must know the adresses as they are created.
	//CompoundState* compoundstates_host;		
	//CompoundNeighborList* compoundneighborlists_host;
	//----------------------------------------------------------------------//
	
	

	// ---------------------------------------------------- Helper functions ---------------------------------------------------- //
	Float3 get3Random() {	// Returns 3 numbers between 0-1
		return Float3(
			rand() % RAND_MAX / (double) RAND_MAX,
			rand() % RAND_MAX / (double)RAND_MAX,
			rand() % RAND_MAX / (double)RAND_MAX
		);
	}

	double random() {
		return rand() % 10000 / 10000.f * 2 - 1.f;
	}
};

