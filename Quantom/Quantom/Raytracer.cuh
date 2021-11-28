#pragma once 


#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <chrono>

#include "QuantomTypes.cuh"
#include "Simulation.cuh"


const double FOCAL_LEN_RATIO = 1;
const double FOCAL_Y_OFFSET = -1;
const int RAYS_PER_DIM = 1000;	// Can't launch kernels if above 1024
const int NUM_RAYS = RAYS_PER_DIM * RAYS_PER_DIM;
const int THREADS_PER_BLOCK = 1024;
const int MAX_RAY_BLOCKS = 80;


class Ray {
public:
	Ray(){}
	Ray(Float3 unit_vector, Float3 origin);
	__device__ void reset();

	__device__ bool hitsParticle(Float3* particle_center, double particle_radius);
	__device__ void searchCompound(CompoundState* compoundstate, Box* box, int i);	// I is for bugfinding.
	__device__ bool searchSolvent(Float3* pos, Box* box, int solvent_index);
	__device__ double distToPoint(Float3 point);
	

	Float3 origin;
	Float3 unit_vector;



	// Changing with each render step
	double closest_collision = 0;
	int atom_type = -1;
	bool log_particle = false;


private:
	__device__ double distToSphereIntersect(Float3* particle_center, double particle_radius);

	
};



class Raytracer {
public:
	Raytracer(){}
	Raytracer(Simulation* simulation, bool verbose);
	uint8_t* render(Simulation* simulation);
	//void renderKernel(Ray* rayptr, uint8_t* image, Float3 focalpoint);

private:
	void initRays();
	Ray* rayptr;
	Float3 focalpoint;


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