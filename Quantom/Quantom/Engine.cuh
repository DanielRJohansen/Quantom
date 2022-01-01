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

	
#include <fstream>	// TEMP


__global__ void forceKernel(Box* box);
__global__ void solventForceKernel(Box* box);
__global__ void initKernel(Box* box);	// For now, just initializes previous potential energy





class Engine {
public:
	Engine(){}

	Simulation* prepSimulation(Simulation* simulation, Compound* main_molecule);
	//double* getDatabuffer();
	bool testFunction();
	void step();
	double* analyzeEnergy();


	Int3 timings = Int3(0, 0, 0);


	// ----- Functions used by analyzer aswell ----- //
	//__device__ Float3 computeLJForces(Box* box, Compound* compound, CompoundNeighborList* neighborlist, 
		//CompoundState* self_state, CompoundState* neighborstate_buffer, Float3* utility_buffer);





	~Engine() {
		//for (int i = 0; i < N_STREAMS; i++)
			//cudaStreamDestroy(stream[i]);
	}
private:
	BoxBuilder boxbuilder;
	Simulation* simulation;


	// HOST FUNCTIONS //
	void updateNeighborLists();






	// ################################# VARIABLES AND ARRAYS ################################# //

	int testval = 0;


	// Simulation variables
	//cudaStream_t stream[N_STREAMS];
	dim3 gridblock_size;
	int threads_per_gridblock;

	bool finished = false;
	int sim_blocks;

	double block_dist;
	int bpd;
	double box_base;				// Of box (left, back, down-most), is negative!
	double block_center_base;	// Including that edge blocks focus area is halfway outside the box
	cudaError_t cuda_status;

};


// FUNCTIONS IN THIS CLASS CAN BE CALLED FROM OTHER FILES
/*
class DeviceFunctionCollection {
public:
	__device__ __host__ DeviceFunctionCollection(){}

	__device__ static double getAngle(Float3 v1, Float3 v2);
	__device__ static void determineMoleculerHyperposOffset(Float3* utility_buffer, Float3* compund_center, 
		Float3* neighbor_center);	

	__device__ static Float3 calcLJForce(Float3* pos0, Float3* pos1);

	__device__ static Float3 calcPairbondForce(Float3* self_pos, Float3* other_pos, double* reference_dist);
	__device__ static Float3 calcAngleForce(CompoundState* statebuffer, AngleBond* anglebond);

	__device__ static Float3 computeLJForces(Box* box, Compound* compound, CompoundNeighborList* neighborlist,
		CompoundState* self_state, CompoundState* neighborstate_buffer, Float3* utility_buffer);
	__device__ static Float3 computePairbondForces(Compound* compound, CompoundState* self_state);
	__device__ static Float3 computeAnglebondForces(Compound* compound, CompoundState* self_state);

};
*/