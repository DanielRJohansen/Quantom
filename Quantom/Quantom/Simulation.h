#pragma once

//#include "QuantomTypes.cuh"
#include "Bodies.cuh"

class Simulation {


public:
	Double3 box_size = Double3(10, 10, 10);
	int n_steps = 1000000;


	
	int n_bodies = 500;
	SimBody* bodies;


};