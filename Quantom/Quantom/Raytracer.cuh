#pragma once 


#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "QuantomTypes.cuh"
#include "Simulation.cuh"


const int FOCAL_LEN_RATIO = 10;
const int RAYS_PER_DIM = 1000;	// Can't launch kernels if above 1024
const int NUM_RAYS = RAYS_PER_DIM * RAYS_PER_DIM;
const int THREADS_PER_BLOCK = 1024;
const int MAX_RAY_BLOCKS = 20;


class Ray {
public:
	Ray(){}
	Ray(Double3 unit_vector);

	__device__ void findBlockHits(Box* box, Double3 focalpoint);
	__device__ bool hitsBody(SimBody* body, Double3 focalpoint);


	Double3 unit_vector;

	int block_indexes[MAX_RAY_BLOCKS];	// BAD. each ray can only hit 20 boxes, before it fades



private:


	__device__ bool hitsBlock(Block* Block, Double3 focalpoint);
};



class Raytracer {
public:
	Raytracer(){}
	Raytracer(Simulation* simulation);
	uint8_t* render(Simulation* simulation);
	//void renderKernel(Ray* rayptr, uint8_t* image, Double3 focalpoint);

private:
	void initRays();
	Ray* rayptr;
	Double3 focalpoint;

	cudaError_t cuda_status;



	// HELPER FUNCTIONS;
	void setGPU() {
		cuda_status = cudaSetDevice(0);
		if (cuda_status != cudaSuccess) {
			fprintf(stderr, "cudaSetDevice failed!");
			exit(1);
		}
	}
};