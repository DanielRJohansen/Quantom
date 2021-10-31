#pragma once 


#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <chrono>

#include "QuantomTypes.cuh"
#include "Simulation.cuh"


const float FOCAL_LEN_RATIO = 1;
const float FOCAL_Y_OFFSET = -1;
const int RAYS_PER_DIM = 1000;	// Can't launch kernels if above 1024
const int NUM_RAYS = RAYS_PER_DIM * RAYS_PER_DIM;
const int THREADS_PER_BLOCK = 1024;
const int MAX_RAY_BLOCKS = 80;


class Ray {
public:
	Ray(){}
	Ray(Float3 unit_vector, Float3 origin);
	__device__ void reset();

	__device__ bool hitsParticle(CompactParticle* particle, float particle_radius);
	__device__ bool moleculeCollisionHandling(Particle* particle, MoleculeLibrary* mol_library, uint8_t* image);

	__device__ void searchCompound(CompoundState* compoundstate, Box* box);

	__device__ float distToPoint(Float3 point);
	

	Float3 origin;
	Float3 unit_vector;



	// Changing with each render step
	float closest_collision = 0;
	int atom_type = -1;
	


private:
	__device__ float distToSphereIntersect(CompactParticle* particle, float particle_radius);

	
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