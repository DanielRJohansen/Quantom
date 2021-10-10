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

	__device__ void findBlockHits(Box* box, Float3 focalpoint);
	__device__ bool hitsParticle(Particle* particle);
	__device__ bool searchBlock(Block* block, uint8_t* image);
	__device__ bool moleculeCollisionHandling(Particle* particle, MoleculeLibrary* mol_library, uint8_t* image);
	__device__ bool hitsBlock(Float3* blockmin, Float3* blockmax, Float3* focalpoint);


	__host__ __device__ float distToPoint(Float3 point) {
		Float3 far_ray_point = origin + unit_vector * 99999999;	////BAAAAAAAAAAAAAD
		return (
			((far_ray_point - origin).cross(origin - point)).len()
			/
			(far_ray_point - origin).len()
			);
	}

	Float3 origin;
	Float3 unit_vector;

	int block_indexes[MAX_RAY_BLOCKS];	// BAD. each ray can only hit 20 boxes, before it fades



private:
	//__device__ float Ray::distToSphereIntersect(Atom* atom);
	__device__ float Ray::distToSphereIntersect(Particle* particle);

	__device__ bool hitsBlock(Block* Block, Float3 focalpoint);
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