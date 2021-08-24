#pragma once

//#include "QuantomTypes.cuh"
#include "Bodies.cuh"

__global__ class Simulation {


public:
	double box_size = 10;	//nm
	int n_steps = 1000000;


	
	int n_bodies = 500;
	SimBody* bodies;


};