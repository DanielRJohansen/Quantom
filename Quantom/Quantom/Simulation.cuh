#pragma once

//#include "QuantomTypes.cuh"
#include "Bodies.cuh"



constexpr double BOX_LEN = 4.f;	// Multiple of FOCUS_LEN please!

constexpr double BOX_LEN_HALF = BOX_LEN/2.f;
constexpr double BOX_BASE = -BOX_LEN_HALF;

constexpr double FOCUS_LEN = 2;
constexpr double BLOCK_LEN = FOCUS_LEN * 2;	//nm
constexpr double FOCUS_LEN_HALF = BLOCK_LEN / 4.f;



const int LOGBLOCK = 0;
const int LOGTHREAD = 50;
const int LOGTYPE = 0;	// 0 for solvent, 1 for compound
//const int N_BODIES_START = 40;
const int N_SOLVATE_MOLECULES = 5*5*5;// 60;
//const int N_SOLVATE_MOLECULES = 3;


const int PARTICLES_PER_COMPOUND = 3;
const int BLOCKS_PER_SM = 512;
//const int GRIDBLOCKS_PER_BODY = 16;
//const int THREADS_PER_GRIDBLOCK = MAX_BLOCK_BODIES / GRIDBLOCKS_PER_BODY;
const int N_STREAMS = 50;			// 68 total, 0 is general purpose, 1 is for rendering.


constexpr double WARN_FORCE = 80'000;
constexpr double END_SIM_FORCE = 10'500'000;
const int LOG_P_ID = 17;

const int MAX_COMPOUNDS = 0xFF;
const int MAX_SOLVENTS = 0xFFFF;
constexpr double CUTOFF = 5.0f;	//nm/


class Box {
public:

	Compound* compounds;

	uint32_t n_compounds = 0;
	uint16_t n_solvents = 0;
	uint32_t total_particles = 0;
	//uint32_t n_solvents = 0;

	RenderMolecule rendermolecule;	// Not proud, TEMP

	// These are shared for all compounds, MUST be allocated before adding any compounds to box, so not in moveToDevice //
	CompoundState* compound_state_array;	
	CompoundState* compound_state_array_next;

	CompoundNeighborList* compound_neighborlist_array;	// First (MAX_COMPUNDS) lists belonging to compounds, then solvents
	SolventNeighborList* solvent_neighborlist_array;	// First (MAX_COMPUNDS) lists belonging to compounds, then solvents
	//------------------------------------//

	Solvent* solvents;
	Solvent* solvents_next;



	double* outdata;			// Temp, for longging values to whatever
	uint32_t step = 0;
	double dt;


	double* potE_buffer;		// For total energy summation
	Float3* trajectory;


	void moveToDevice() {	// Loses pointer to RAM location!
		compounds = genericMoveToDevice(compounds, n_compounds);

		solvents = genericMoveToDevice(solvents, MAX_SOLVENTS);

		compound_state_array = genericMoveToDevice(compound_state_array, MAX_COMPOUNDS);



		compound_neighborlist_array = genericMoveToDevice(compound_neighborlist_array, MAX_COMPOUNDS + MAX_SOLVENTS);
		solvent_neighborlist_array = genericMoveToDevice(solvent_neighborlist_array, MAX_COMPOUNDS + MAX_SOLVENTS);


		


		cudaDeviceSynchronize();

		printf("Box transferred to device\n\n");
	}
	
};

	


__global__ class Simulation {


public:
	Simulation() {
		box = new Box;
	}

	void moveToDevice() {
		box = genericMoveToDevice(box, 1);
		printf("Simulation ready for device\n");
	}

	bool finished = false;
	//int step = 0;


	double box_size = BOX_LEN;	//nm
	int blocks_per_dim;
	int n_steps = 200;
	//int n_steps = 3000;
	const double dt = 2 * 1e-6;		// ns, so first val corresponds to fs
	int steps_per_render = 50;

	//int n_bodies = N_BODIES_START;
	Box* box;

	~Simulation() {

	}



};







