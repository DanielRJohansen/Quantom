#pragma once

//#include "QuantomTypes.cuh"
#include "Bodies.cuh"



constexpr float BOX_LEN = 6.0;	// Multiple of FOCUS_LEN please!

constexpr float BOX_LEN_HALF = BOX_LEN/2.f;
constexpr float BOX_BASE = -BOX_LEN_HALF;

constexpr float FOCUS_LEN = 2;
constexpr float BLOCK_LEN = FOCUS_LEN * 2;	//nm
constexpr float FOCUS_LEN_HALF = BLOCK_LEN / 4.f;

constexpr float edgeforce_scalar = 100000.f;
constexpr float edgeforce_scalar_half = edgeforce_scalar/2.f;


//constexpr auto CUTOFF_LEN = 0.8f;		// nm
//constexpr float BLOCK_OVERLAP = 0.3f;	// nm, must be > 2* vdw radius of largest atom.

const int MAX_FOCUS_BODIES = 32;
const int MAX_NEAR_BODIES = 256 - MAX_FOCUS_BODIES;
//constexpr float SOLOBLOCK_DIST = BLOCK_LEN - BLOCK_OVERLAP;


const int INDEXA = 100900;
//const int N_BODIES_START = BOX_LEN*BOX_LEN*BOX_LEN/(FOCUS_LEN*FOCUS_LEN*FOCUS_LEN) * 25;
const int N_BODIES_START = 14;
const int BLOCKS_PER_SM = 512;
//const int GRIDBLOCKS_PER_BODY = 16;
//const int THREADS_PER_GRIDBLOCK = MAX_BLOCK_BODIES / GRIDBLOCKS_PER_BODY;
const int N_STREAMS = 50;			// 68 total, 0 is general purpose, 1 is for rendering.




// USEFUL VALUES



struct Block {	// All boxes are cubic

	__host__ __device__ Block() {}
	__host__ __device__ Block(Float3 center) : center(center) {}
	__host__ __device__ bool isInBLock(Float3 point);

	__host__ bool addParticle(Particle* particle) {			// ONLY USED FOR INITIATION 
		if (n_particles == MAX_FOCUS_BODIES-8) {
			//printf("Too many particles for this block!");
			return false;
		}
		//printf("Added\n");
		focus_particles[n_particles] = *particle;
		n_particles++;
		return true;
	}


	Float3 center;

	Particle focus_particles[MAX_FOCUS_BODIES];
	Particle near_particles[MAX_NEAR_BODIES];

	int n_particles = 0;		//  Only used when loading the block
	bool edge_block = false;

	//signed char edge_type[3] = { 0,0,0 }; // -1 negative edge, 0 non-edge, 1 pos-edge, xyz

};

struct AccessPoint {
	Particle particles[MAX_FOCUS_BODIES];
};

class Box {	// Should each GPU block have a copy of the box?
public:
	int n_blocks;
	Block* blocks;
	
	AccessPoint* accesspoint;

	int blocks_per_dim;

	Particle* particles;
	uint32_t n_particles = 0;
	//CompactParticle* all_particles;	// Used for communation between inter and intra-molecular kernels
	
	Compound_H2O* compounds;
	uint32_t n_compounds = 0;

	uint32_t n_pairbonds = 0;	//Need to record this so we can avoid LJ pot for bonded particles

	void moveToDevice() {	// Loses pointer to RAM location!
		//printf("Block 38: %.1f %.1f %.1f\n", blocks[38].center.x, blocks[38].center.y, blocks[38].center.z);

		Block* blocks_temp;
		int blocks_bytesize = n_blocks * sizeof(Block);
		cudaMallocManaged(&blocks_temp, blocks_bytesize);
		cudaMemcpy(blocks_temp, blocks, blocks_bytesize, cudaMemcpyHostToDevice);
		cudaDeviceSynchronize();
		delete blocks;
		blocks = blocks_temp;


		Particle* particles_temp;
		int bytesize = n_particles * sizeof(Particle);
		cudaMalloc(&particles_temp, bytesize);
		cudaMemcpy(particles_temp, particles, bytesize, cudaMemcpyHostToDevice);
		delete particles;
		particles = particles_temp;

		Compound_H2O* compounds_temp;
		bytesize = n_compounds * sizeof(Compound_H2O);
		cudaMalloc(&compounds_temp, bytesize);
		cudaMemcpy(compounds_temp, compounds, bytesize, cudaMemcpyHostToDevice);
		delete compounds;
		compounds = compounds_temp;


		AccessPoint* ap_temp = new AccessPoint[n_blocks];
		for (int i = 0; i < n_blocks; i++) {
			ap_temp[i] = AccessPoint();
		}
		int ap_bytesize = n_blocks * sizeof(AccessPoint);
		cudaMalloc(&accesspoint, ap_bytesize);
		cudaMemcpy(accesspoint, ap_temp, ap_bytesize, cudaMemcpyHostToDevice);
		delete ap_temp;
		//printf("Block 38: %.1f %.1f %.1f\n", blocks[38].center.x, blocks[38].center.y, blocks[38].center.z);
		//printf("Block 38: %.1f %.1f %.1f\n", blocks_temp[38].center.x, blocks_temp[38].center.y, blocks_temp[38].center.z);
		printf("Box transferred to device\n\n");
	}



};

	


__global__ class Simulation {


public:
	Simulation() {}
	Simulation(MoleculeLibrary* mol_library) : mol_library(mol_library) {
		box = new Box;
	}

	void moveToDevice() {
		box->moveToDevice();

		Box* box_temp;
		cudaMallocManaged(&box_temp, sizeof(Box));
		cudaMemcpy(box_temp, box, sizeof(Box), cudaMemcpyHostToDevice);
		delete box;
		box = box_temp;

		MoleculeLibrary* lib_temp;
		cudaMallocManaged(&lib_temp, sizeof(MoleculeLibrary));
		cudaMemcpy(lib_temp, mol_library, sizeof(MoleculeLibrary), cudaMemcpyHostToDevice);
		delete mol_library;
		mol_library = lib_temp;

		
		printf("Simulation ready for device\n");
	}

	bool finished = false;
	int step = 0;


	float box_size = BOX_LEN;	//nm
	int blocks_per_dim;
	int n_steps = 60000;

	const double dt = 1 *	10.0e-6;		// ns, so first val corresponds to fs
	int steps_per_render = 50;

	int n_bodies = N_BODIES_START;
	Box* box;
	//SimBody* bodies;	// The bodies of each block is only total copy, not a pointer to its corresponding body here!
	MoleculeLibrary* mol_library;

	~Simulation() {
		//delete box->blocks;
		//delete box;
	}



};







