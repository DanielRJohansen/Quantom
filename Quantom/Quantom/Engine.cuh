#pragma once

#include <iostream>
#include <chrono>

#include "Bodies.cuh"
#include "Simulation.cuh"

#include <cuda.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <device_functions.h>
#include <cuda_runtime_api.h>

	

__global__ void stepKernel(Simulation* simulation, int offset);
__global__ void updateKernel(Simulation* simulation, int offset);


class Engine {
public:
	Engine(){}

	Simulation* prepSimulation(Simulation* simulation);

	int countBodies();
	void step();


	Int3 timings = Int3(0, 0, 0);


	~Engine() {
		for (int i = 0; i < N_STREAMS; i++)
			cudaStreamDestroy(stream[i]);
	}
private:
	Simulation* simulation;
	int fillBox();		// Returns # of bodies placed

	int initBlocks();	// returns # of blocks created
	void linkBlocks();
	void prepareEdgeBlocks();

	void placeBody(SimBody* body);
	void prepareCudaScheduler();

	Simulation* simToDevice();





	// Helper functions
	/*Int3 posToBlockIndex(Float3* pos) {
		return Int3(
			floor((pos->x + simulation->box_size / 2.f) / simulation->box_size * simulation->blocks_per_dim),
			floor((pos->y + simulation->box_size / 2.f) / simulation->box_size * simulation->blocks_per_dim),
			floor((pos->z + simulation->box_size / 2.f) / simulation->box_size * simulation->blocks_per_dim)
			);
	}*/
	Int3 posToBlockIndex(Float3* pos) {
		return Int3(
			floor((pos->x - box_base) / (bpd * block_dist) * simulation->blocks_per_dim),
			floor((pos->y - box_base) / (bpd * block_dist) * simulation->blocks_per_dim),
			floor((pos->z - box_base) / (bpd * block_dist) * simulation->blocks_per_dim)
		);
	}
	inline int block3dIndexTo1dIndex(Int3 index_3d) {// NOTICE THAT X IS THE "GRANDPARENT" ITERATOR
		return (index_3d.x + 
			index_3d.y * bpd + 
			index_3d.z * bpd*bpd);
	}






	// Simulation variables
	cudaStream_t stream[N_STREAMS];
	dim3 gridblock_size;
	int threads_per_gridblock;

	bool finished = false;
	int sim_blocks;

	float block_dist;
	int bpd;
	float box_base;				// Of box (left, back, down-most), is negative!
	float block_center_base;	// Including that edge blocks focus area is halfway outside the box
	cudaError_t cuda_status;

};