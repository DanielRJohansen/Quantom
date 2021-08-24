#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "math.h"



__global__ struct Int3 {
	Int3() {}
	Int3(int x, int y, int z) : x(x), y(y), z(z) {}


	int x, y, z;
};

__host__ __device__ struct Double3 {
	__host__ __device__ Double3() {}
	__host__ __device__ Double3(double x, double y, double z) : x(x), y(y), z(z) {}

	__host__ __device__ inline Double3 operator * (const double a) const { return Double3(x * a, y * a, z * a); }
	__host__ __device__ inline Double3 operator + (const Double3 a) const { return Double3(x + a.x, y + a.y, z + a.z); }
	__host__ __device__ inline Double3 operator - (const Double3 a) const { return Double3(x - a.x, y - a.y, z - a.z); }


	__host__ __device__ inline double length();




	//Utilities
	//void print()

	double x, y, z;
};








__device__ class Block {	// All boxes are cubic
public:
	Block() {}

	bool isInBLock(Double3 point);

	Double3 center;
	double length;
private:

};

__global__ class Box {
public:
	int n_blocks;
	Block* blocks;

	void moveToDevice();	// Loses pointer to RAM location!
};