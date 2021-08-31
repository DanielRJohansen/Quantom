#pragma once

//#include "QuantomTypes.cuh"
#include "Bodies.cuh"


const int MAX_BLOCK_BODIES = 512;
//const double BLOCK_LEN = 10; //nm
constexpr auto BLOCK_LEN_CUDA = 1.0; //nm
constexpr auto BOX_LEN_CUDA = 10.0;
//#define BLOCK_LEN 1.0;

const int INDEXA = 999;
const int N_BODIES_START = 10000;






class Block {	// All boxes are cubic
public:
	__host__ __device__ Block() {}
	__host__ __device__ Block(Double3 center) : center(center) {}
	__host__ __device__ bool isInBLock(Double3 point);

	Double3 center;
	SimBody bodies[MAX_BLOCK_BODIES];
	int n_bodies = 0;

private:

};

class Box {
public:
	int n_blocks;
	Block* blocks;

	void moveToDevice() {	// Loses pointer to RAM location!
		//printf("Block 38: %.1f %.1f %.1f\n", blocks[38].center.x, blocks[38].center.y, blocks[38].center.z);

		Block* blocks_temp;
		int blocks_bytesize = n_blocks * sizeof(Block);
		cudaMallocManaged(&blocks_temp, blocks_bytesize);
		cudaMemcpy(blocks_temp, blocks, blocks_bytesize, cudaMemcpyHostToDevice);
		cudaDeviceSynchronize();
		delete blocks;
		blocks = blocks_temp;
		//printf("Block 38: %.1f %.1f %.1f\n", blocks[38].center.x, blocks[38].center.y, blocks[38].center.z);
		//printf("Block 38: %.1f %.1f %.1f\n", blocks_temp[38].center.x, blocks_temp[38].center.y, blocks_temp[38].center.z);

	}
};







__global__ class Simulation {


public:
	Simulation() {
		box = new Box;
	}
	void moveToDevice() {
		box->moveToDevice();

		Box* box_temp;
		cudaMallocManaged(&box_temp, sizeof(Box));
		cudaMemcpy(box_temp, box, sizeof(Box), cudaMemcpyHostToDevice);
		delete box;
		box = box_temp;
		
		
		printf("Simulation ready for device\n");
	}


	double box_size = BOX_LEN_CUDA;	//nm
	int blocks_per_dim;
	int n_steps = 1000000;


	
	int n_bodies = N_BODIES_START;
	Box* box;
	SimBody* bodies;	// The bodies of each block is only total copy, not a pointer to its corresponding body here!

	~Simulation() {
		delete box->blocks;
		delete box;
	}



};







