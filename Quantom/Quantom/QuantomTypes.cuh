#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"




__global__ struct Int3 {
	Int3() {}
	Int3(int x, int y, int z) : x(x), y(y), z(z) {}


	int x, y, z;
};

__global__ struct Double3 {
	Double3() {}
	Double3(double x, double y, double z) : x(x), y(y), z(z) {}

	inline Double3 operator * (const double a) const { return Double3(x * a, y * a, z * a); }
	inline Double3 operator + (const Double3 a) const { return Double3(x + a.x, y + a.y, z + a.z); }
	inline Double3 operator - (const Double3 a) const { return Double3(x - a.x, y - a.y, z - a.z); }







	//Utilities
	//void print()

	double x, y, z;
};
