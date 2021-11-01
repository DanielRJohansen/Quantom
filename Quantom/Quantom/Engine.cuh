#pragma once

#include <iostream>
#include <chrono>

#include "Bodies.cuh"
#include "Simulation.cuh"
#include "BoxBuilder.cuh"

#include <cuda.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <device_functions.h>
#include <cuda_runtime_api.h>

	


__global__ void forceKernel(Box* box);


class Engine {
public:
	Engine(){}

	Simulation* prepSimulation(Simulation* simulation);

	void step();

	Int3 timings = Int3(0, 0, 0);


	~Engine() {
		for (int i = 0; i < N_STREAMS; i++)
			cudaStreamDestroy(stream[i]);
	}
private:
	BoxBuilder boxbuilder;
	Simulation* simulation;

	void updateNeighborLists();






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