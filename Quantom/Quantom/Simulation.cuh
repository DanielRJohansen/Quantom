#pragma once

//#include "QuantomTypes.cuh"
#include "Bodies.cuh"


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


	double box_size = 10;	//nm
	double blocks_per_dim;
	int n_steps = 1000000;


	
	int n_bodies = 500;
	Box* box;
	SimBody* bodies;

	~Simulation() {
		delete box->blocks;
		delete box;
	}



};