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

	
#include <thread>

#include <fstream>	// TEMP


__global__ void forceKernel(Box* box);
__global__ void solventForceKernel(Box* box);
__global__ void initKernel(Box* box);	// For now, just initializes previous potential energy



struct NListDataCollection {
	//NListDataCollection() {}
	NListDataCollection(Simulation* simulation) {
		n_compounds = simulation->n_compounds;
		n_solvents = simulation->n_solvents;
		compoundstates = new CompoundState[n_compounds];
		solvents = new Solvent[simulation->n_solvents];
		compound_neighborlists = new NeighborList[MAX_COMPOUNDS];
		solvent_neighborlists = new NeighborList[MAX_SOLVENTS];
		cudaMemcpy(compound_neighborlists, simulation->box->compound_neighborlists, sizeof(NeighborList) * n_compounds, cudaMemcpyDeviceToHost);
		cudaMemcpy(solvent_neighborlists, simulation->box->solvent_neighborlists, sizeof(NeighborList) * n_solvents, cudaMemcpyDeviceToHost);
	}
	void compressPositionData() {
		for (int i = 0; i < n_compounds; i++) {
			compound_key_positions[i] = compoundstates[i].positions[0];
		}
		for (int i = 0; i < n_compounds; i++) {
			solvent_positions[i] = solvents[i].pos;
		}
	}
	int n_compounds;
	int n_solvents;

	// I guess this is not critical but temp, needed to load pos device->host
	CompoundState* compoundstates;
	Solvent* solvents;

	Float3 compound_key_positions[MAX_COMPOUNDS];
	Float3 solvent_positions[MAX_SOLVENTS];

	// These are loaded before simulaiton start. Kept on host, and copied to device each update.
	NeighborList* compound_neighborlists;
	NeighborList* solvent_neighborlists;
};



class Engine {
public:
	Engine();
	Engine(Simulation* simulation);


	void deviceMaster();
	void hostMaster();

	//double* getDatabuffer();
	bool testFunction();
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
	Simulation* simulation;


	// -------------------------------------- GPU LOAD -------------------------------------- //
	void step();

	// -------------------------------------- CPU LOAD -------------------------------------- //
	void offloadPositionData(Simulation* simulation);
	void onloadNeighborlists();
	static void updateNeighborLists(Simulation* simulation, NListDataCollection* nlist_data_collection, 
		volatile bool* finished, int* timing);	// thread worker, can't own engine object, thus pass ref
	static void cullDistantNeighbors(NListDataCollection* nlist_data_collection);
	NListDataCollection* nlist_data_collection;

	void offloadLoggingData();



	int prev_nlist_update_step = 0;
	volatile bool updated_neighborlists_ready = 0;

	// -------------------------------------- HELPERS -------------------------------------- //

	void genericErrorCheck(const char* text) {
		cudaError_t cuda_status = cudaGetLastError();
		if (cuda_status != cudaSuccess) {
			fprintf(stderr, text);
			exit(1);
		}
	}




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