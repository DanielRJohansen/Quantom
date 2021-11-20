#pragma once

//#include "QuantomTypes.cuh"
#include "Bodies.cuh"



constexpr float BOX_LEN = 2.0;	// Multiple of FOCUS_LEN please!

constexpr float BOX_LEN_HALF = BOX_LEN/2.f;
constexpr float BOX_BASE = -BOX_LEN_HALF;

constexpr float FOCUS_LEN = 2;
constexpr float BLOCK_LEN = FOCUS_LEN * 2;	//nm
constexpr float FOCUS_LEN_HALF = BLOCK_LEN / 4.f;

constexpr float edgeforce_scalar = 100000.f;
constexpr float edgeforce_scalar_half = edgeforce_scalar/2.f;


//constexpr auto CUTOFF_LEN = 0.8f;		// nm
//constexpr float BLOCK_OVERLAP = 0.3f;	// nm, must be > 2* vdw radius of largest atom.

const int MAX_FOCUS_BODIES = 120;
const int MAX_NEAR_BODIES = 256 - MAX_FOCUS_BODIES;
//constexpr float SOLOBLOCK_DIST = BLOCK_LEN - BLOCK_OVERLAP;


const int LOGBLOCK = 0;
const int LOGTHREAD = 0;

//const int N_BODIES_START = 40;
const int N_SOLVATE_MOLECULES = 30;// 60;
const int BLOCKS_PER_SM = 512;
//const int GRIDBLOCKS_PER_BODY = 16;
//const int THREADS_PER_GRIDBLOCK = MAX_BLOCK_BODIES / GRIDBLOCKS_PER_BODY;
const int N_STREAMS = 50;			// 68 total, 0 is general purpose, 1 is for rendering.


constexpr float WARN_FORCE = 80'000;
constexpr float END_SIM_FORCE = 10'500'000;
const int LOG_P_ID = 17;



class Box {
public:

	Compound* compounds;

	uint32_t n_compounds = 0;
	//uint32_t n_solvents = 0;

	RenderMolecule rendermolecule;	// Not proud, TEMP

	// These are shared for all compounds, MUST be allocated before adding any compounds to box, so not in moveToDevice //
	CompoundState* compound_state_array;	
	CompoundState* compound_state_array_next;
	CompoundNeighborList* compound_neighborlist_array;
	//------------------------------------//

	float* outdata;	// Temp, for longging values to whatever
	uint32_t step = 0;
	float dt;


	float* data_buffer;		// also temp, for total energy summation

	Float3* trajectory;


	void moveToDevice() {	// Loses pointer to RAM location!
		Compound* compounds_temp;
		int bytesize = n_compounds * sizeof(Compound);
		cudaMallocManaged(&compounds_temp, bytesize);
		cudaMemcpy(compounds_temp, compounds, bytesize, cudaMemcpyHostToDevice);
		delete compounds;
		compounds = compounds_temp;



		cudaMallocManaged(&outdata, sizeof(float) * 10 * 10000);	// 10 data streams for 10k steps. 1 step at a time.


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
		Box* box_temp;
		cudaMallocManaged(&box_temp, sizeof(Box));
		cudaMemcpy(box_temp, box, sizeof(Box), cudaMemcpyHostToDevice);
		delete box;
		box = box_temp;

		printf("Simulation ready for device\n");
	}

	bool finished = false;
	int step = 0;


	float box_size = BOX_LEN;	//nm
	int blocks_per_dim;
	int n_steps = 4000;

	const float dt = 1 * 1e-6;		// ns, so first val corresponds to fs
	int steps_per_render = 40;

	//int n_bodies = N_BODIES_START;
	Box* box;

	~Simulation() {

	}



};







