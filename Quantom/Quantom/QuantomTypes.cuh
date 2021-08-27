#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "math.h"

const int BLOCK_LEN = 1; //nm





struct Int3 {
	__host__ __device__ Int3() {}
	__host__ __device__ Int3(int x, int y, int z) : x(x), y(y), z(z) {}


	int x=0, y=0, z = 0;
};

struct Double3 {
	__host__ __device__ Double3() {}
	__host__ __device__ Double3(double x, double y, double z) : x(x), y(y), z(z) {}

	__host__ __device__ inline Double3 operator * (const double a) const { return Double3(x * a, y * a, z * a); }
	__host__ __device__ inline Double3 operator + (const Double3 a) const { return Double3(x + a.x, y + a.y, z + a.z); }
	__host__ __device__ inline Double3 operator - (const Double3 a) const { return Double3(x - a.x, y - a.y, z - a.z); }


	__host__ __device__ inline double len() {return (double)sqrtf(x * x + y * y + z * z); }




	//Utilities
	//void print()

	double x, y, z;
};








class Block {	// All boxes are cubic
public:
	__host__ __device__ Block() {}
	__host__ __device__ Block(Double3 center) : center(center){}
	__host__ __device__ bool isInBLock(Double3 point);

	Double3 center;

private:

};

class Box {
public:
	int n_blocks;
	Block* blocks;

	void moveToDevice() {	// Loses pointer to RAM location!
		Block* blocks_temp;
		int blocks_bytesize = n_blocks * sizeof(Block);
		cudaMallocManaged(&blocks_temp, blocks_bytesize);
		cudaMemcpy(blocks_temp, blocks, sizeof(Box), cudaMemcpyHostToDevice);
		delete blocks;
		blocks = blocks_temp;
	}	
};