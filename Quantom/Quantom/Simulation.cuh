#pragma once

//#include "QuantomTypes.cuh"
#include "Bodies.cuh"



constexpr double BOX_LEN = 6.f;		// Must be > twice the len of largest compound
constexpr double BOX_LEN_HALF = BOX_LEN/2.f;



const int STEPS_PER_NLIST_UPDATE = 50;
const int STEPS_PER_LOGTRANSFER = 100;
//const int STEPS_PER_TRAJTRANSFER = 100;
const int STEPS_PER_THERMOSTAT = 100;
const int STEPS_PER_TRAINDATATRANSFER = 40;



const int STEPS_PER_RENDER = 50;


const bool APPLY_THERMOSTAT = true;

#ifndef __linux__
const string MOL_FOLDER = "C:\\PROJECTS\\Quantom\\t4lys\\";
const string OUT_DIR = "C:\\PROJECTS\\Quantom\\LIMANET\\sim_out";
#else
const string MOL_FOLDER = "../Compounds/t4lys/";
const string OUT_DIR = "/home/lima/Desktop/LIMA";
#endif

const int N_SOLVATE_MOLECULES = 800;	// Must not be above 255, as that can't be represented as uint8_t5





constexpr double WARN_FORCE = 80'000;
constexpr double END_SIM_FORCE = 10'500'000;
const int LOG_P_ID = 17;

const int MAX_COMPOUNDS = 0xFF;
const int MAX_SOLVENTS = 0xFFFF;
constexpr float CUTOFF = 0.8f;	//nm/



const int THREADS_PER_SOLVENTBLOCK = 256;	// Must be >= N_SOLVATE_MOLECULES
const int BLOCKS_PER_SOLVENTKERNEL = ceil((float)N_SOLVATE_MOLECULES/(float)THREADS_PER_SOLVENTBLOCK);


const int THREADS_PER_COMPOUNDBLOCK = 128; // Must be >= max comp particles

const int N_LIPID_COPIES = 32;


const int SIMULATION_STEPS = 20000;

//constexpr double SOLVENT_MASS = 18.01528f * 1e-3;	// kg/mol
constexpr double SOLVENT_MASS = 12.0107 * 1e-3;	// kg/mol
constexpr double COMPOUNDPARTICLE_MASS = 12.0107 * 1e-3;

// This goes on Device
class Box {
public:

	Compound* compounds;

	uint32_t n_compounds = 0;
	uint16_t n_solvents = 0;
	uint32_t total_particles_upperbound = 0;

	RenderMolecule rendermolecule;	// Not proud, TEMP

	// These are shared for all compounds, MUST be allocated before adding any compounds to box, so not in moveToDevice //
	CompoundState* compound_state_array;	
	CompoundState* compound_state_array_next;

	NeighborList* compound_neighborlists;
	NeighborList* solvent_neighborlists;
	//------------------------------------//

	Solvent* solvents;
	Solvent* solvents_next;



	uint32_t step = 0;
	double dt;
	bool critical_error_encountered = 0;

	double* potE_buffer;		// For total energy summation
	Float3* traj_buffer;

	double* outdata;			// Temp, for longging values to whatever
	Float3* data_GAN;			// Only works if theres 1 compounds right now.


	float thermostat_scalar = 1.f;


	void moveToDevice() {	// Loses pointer to RAM location!

		int bytes_total = sizeof(Compound) * n_compounds
			+ sizeof(Solvent) * MAX_SOLVENTS * 2
			+ sizeof(CompoundState) * MAX_COMPOUNDS * 2
			+ sizeof(NeighborList) * (MAX_SOLVENTS + MAX_COMPOUNDS);
		printf("Moving %.2f GB to device\n", (float)bytes_total * 1e-9);

		compounds = genericMoveToDevice(compounds, n_compounds);
		solvents = genericMoveToDevice(solvents, MAX_SOLVENTS);
		compound_state_array = genericMoveToDevice(compound_state_array, MAX_COMPOUNDS);


		compound_neighborlists = genericMoveToDevice(compound_neighborlists, MAX_COMPOUNDS);
		solvent_neighborlists = genericMoveToDevice(solvent_neighborlists, MAX_SOLVENTS);


		cudaDeviceSynchronize();

		printf("Box transferred to device\n\n");
	}
	
};

	

// This stays on host
class Simulation {
public:
	Simulation() {
		box = new Box;
	}

	__host__ void moveToDevice() {
		box = genericMoveToDevice(box, 1);
		cudaDeviceSynchronize();
		if (cudaGetLastError() != cudaSuccess) {
			fprintf(stderr, "Error during Simulation Host->Device transfer\n");
			exit(1);
		}
		printf("Simulation ready for device\n");
	}
	__host__ void copyBoxVariables() {
		n_compounds = box->n_compounds;
		n_solvents = box->n_solvents;
		compounds_host = new Compound[n_compounds];
		for (int i = 0; i < n_compounds; i++)
			compounds_host[i] = box->compounds[i];
	}
	__host__ void incStep() {
		step++;
		box->step++;
	}
	__host__ inline int getStep() {
		return step;
	}
	~Simulation() {}

	bool finished = false;
	//int step = 0;


	double* potE_buffer;	// Not really a buffer yet, just one large array that holds full simulation data
	Float3* traj_buffer;
	float* temperature_buffer;
	Float3* traindata_buffer;


	uint32_t total_particles_upperbound = 0;

	
	int n_steps = SIMULATION_STEPS;

	const double dt = 1 * 1e-6;		// ns, so first val corresponds to fs
	int steps_per_render = STEPS_PER_RENDER;
	//int n_bodies = N_BODIES_START;
	Box* box;


	Compound* compounds_host;				// For reading static data, for example during nlist-search

	// Box variable copies, here for ease of access.
	int n_compounds = 0;
	int n_solvents = 0;

	string out_dir = OUT_DIR;
	

private:
	int step = 0;


};

