#pragma once

//#include "QuantomTypes.cuh"
#include "Bodies.cuh"



constexpr float BOX_LEN = 24.0;	// Multiple of 2 please!
constexpr float BLOCK_LEN = 4.0;	//nm
constexpr float FOCUS_LEN = BLOCK_LEN / 2.f;
constexpr float FOCUS_LEN_HALF = BLOCK_LEN / 4.f;




//constexpr auto CUTOFF_LEN = 0.8f;		// nm
//constexpr float BLOCK_OVERLAP = 0.3f;	// nm, must be > 2* vdw radius of largest atom.

const int MAX_FOCUS_BODIES = 16;
const int MAX_NEAR_BODIES = 64 - MAX_FOCUS_BODIES;
//constexpr float SOLOBLOCK_DIST = BLOCK_LEN - BLOCK_OVERLAP;


const int INDEXA = 100900;
const int N_BODIES_START = 800;

const int BLOCKS_PER_SM = 16;
//const int GRIDBLOCKS_PER_BODY = 16;
//const int THREADS_PER_GRIDBLOCK = MAX_BLOCK_BODIES / GRIDBLOCKS_PER_BODY;
const int N_STREAMS = 60;			// 68 total, 0 is general purpose, 1 is for rendering.





struct Block {	// All boxes are cubic

	__host__ __device__ Block() {}
	__host__ __device__ Block(Float3 center) : center(center) {}
	__host__ __device__ bool isInBLock(Float3 point);

	__host__ bool addBody(SimBody* body) {			// ONLY USED FOR INITIATION 
		if (n_bodies == MAX_FOCUS_BODIES) {
			printf("Too many bodies for this block!");
			exit(1);
		}

		focus_bodies[n_bodies] = *body;
		n_bodies++;
		return true;
	}


	Float3 center;

	SimBody focus_bodies[MAX_FOCUS_BODIES];
	SimBody near_bodies[MAX_NEAR_BODIES];

	int n_bodies = 0;		//  Only used when loading the block
	bool edge_block = false;

};

struct AccessPoint {
	SimBody bodies[MAX_FOCUS_BODIES];
};

class Box {	// Should each GPU block have a copy of the box?
public:
	int n_blocks;
	Block* blocks;
	
	AccessPoint* accesspoint;

	int blocks_per_dim;

	void finalizeBlock() {

	}
	void moveToDevice() {	// Loses pointer to RAM location!
		//printf("Block 38: %.1f %.1f %.1f\n", blocks[38].center.x, blocks[38].center.y, blocks[38].center.z);

		Block* blocks_temp;
		int blocks_bytesize = n_blocks * sizeof(Block);
		cudaMallocManaged(&blocks_temp, blocks_bytesize);
		cudaMemcpy(blocks_temp, blocks, blocks_bytesize, cudaMemcpyHostToDevice);
		cudaDeviceSynchronize();
		delete blocks;
		blocks = blocks_temp;


		cudaMalloc(&accesspoint, n_blocks * sizeof(AccessPoint));

		//printf("Block 38: %.1f %.1f %.1f\n", blocks[38].center.x, blocks[38].center.y, blocks[38].center.z);
		//printf("Block 38: %.1f %.1f %.1f\n", blocks_temp[38].center.x, blocks_temp[38].center.y, blocks_temp[38].center.z);
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


	float box_size = BOX_LEN;	//nm
	int blocks_per_dim;
	int n_steps = 1000000;

	const float dt = 0.005;
	
	int n_bodies = N_BODIES_START;
	Box* box;
	SimBody* bodies;	// The bodies of each block is only total copy, not a pointer to its corresponding body here!
	MoleculeLibrary* mol_library;

	~Simulation() {
		//delete box->blocks;
		//delete box;
	}



};







