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

	

__global__ void stepKernel(Simulation* simulation);



class Engine {
public:
	Engine(){}

	Simulation* prepSimulation(Simulation* simulation);


	void step();


	~Engine() {
		for (int i = 0; i < N_STREAMS; i++)
			cudaStreamDestroy(stream[i]);
	}
private:
	Simulation* simulation;
	int fillBox();		// Returns # of bodies placed

	int initBlocks();	// returns # of blocks created
	void linkBlocks();

	void placeBody(SimBody* body);
	void prepareCudaScheduler();

	Simulation* simToDevice();





	// Helper functions
	Int3 posToBlockIndex(Float3* pos) {
		return Int3(
			floor((pos->x + simulation->box_size / 2.f) / simulation->box_size * simulation->blocks_per_dim),
			floor((pos->y + simulation->box_size / 2.f) / simulation->box_size * simulation->blocks_per_dim),
			floor((pos->z + simulation->box_size / 2.f) / simulation->box_size * simulation->blocks_per_dim)
			);
	}
	inline int block3dIndexTo1dIndex(Int3 index_3d) {// NOTICE THAT X IS THE "GRANDPARENT" ITERATOR
		return (index_3d.z + 
			index_3d.y * simulation->blocks_per_dim + 
			index_3d.x * simulation->blocks_per_dim * simulation->blocks_per_dim);
	}






	// Simulation variables
	cudaStream_t stream[N_STREAMS];
	dim3 gridblock_size;
	int threads_per_gridblock;

	bool finished = false;
	int sim_blocks;


	cudaError_t cuda_status;

};