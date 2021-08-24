#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "Simulation.cuh"


const int FOCAL_LEN_RATIO = 10;
const int RAYS_PER_DIM = 1000;
const int NUM_RAYS = RAYS_PER_DIM * RAYS_PER_DIM;
const int THREADS_PER_BLOCK = 1024;
const int MAX_RAY_BLOCKS = 20;


__global__ class Ray {
public:
	Ray(){}
	Ray(Double3 unit_vector);

	__device__ void findBlockHits(Box* box, Double3 focalpoint);



private:
	Double3 unit_vector;

	int block_indexes[MAX_RAY_BLOCKS];	// BAD. each ray can only hit 20 boxes, before it fades

	__device__ bool hitsBlock(Block* Block, Double3 focalpoint);

};



__global__ class Raytracer {
public:
	Raytracer(){}
	Raytracer(Simulation* simulation);


private:
	void initRays();
	Ray* rayptr;
};