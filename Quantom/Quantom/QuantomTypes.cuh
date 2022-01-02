#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <cuda.h>
#include <device_atomic_functions.h>
#include <device_functions.h>

#include "math.h"
#include <iostream>

#include <string>
#include <fstream>

constexpr double PI = 3.14159;



struct Int3 {
	__host__ __device__ Int3() {}
	__host__ __device__ Int3(int x, int y, int z) : x(x), y(y), z(z) {}


	__host__ __device__ inline Int3 operator + (const Int3 a) const { return Int3(x + a.x, y + a.y, z + a.z); }
	__host__ __device__ inline Int3 operator - (const Int3 a) const { return Int3(x - a.x, y - a.y, z - a.z); }


	int x=0, y=0, z = 0;
};

struct Float3 {
	__host__ __device__ Float3() {}
	__host__ __device__ Float3(double x, double y, double z) : x(x), y(y), z(z) {}
	__host__ __device__ Float3(double* a) { x = a[0]; y = a[1]; z = a[2]; }

	__host__ __device__ inline Float3 operator * (const double a) const { return Float3(x * a, y * a, z * a); }
	__host__ __device__ inline Float3 operator * (const Float3 a) const { return Float3(x * a.x, y * a.y, z * a.z); }
	__host__ __device__ inline Float3 operator + (const Float3 a) const { return Float3(x + a.x, y + a.y, z + a.z); }
	__host__ __device__ inline Float3 operator - (const Float3 a) const { return Float3(x - a.x, y - a.y, z - a.z); }
	__host__ __device__ inline bool operator == (const Float3 a) const { return (a.x == x && a.y == y && a.z == z); }
	__host__ __device__ inline void operator += (const Float3 a) { x += a.x; y += a.y; z += a.z; }
	__host__ __device__ inline void operator *= (const double a) { x *= a; y *= a; z *= a; }

	__host__ __device__ Float3 norm() {
		double l = len();
		if (l)
			return *this * (1.f / l); 
		return Float3(0, 0, 0);
	}
	__host__ __device__ Float3 square() {return Float3(x * x, y * y, z * z);}
	__host__ __device__ inline double len() {return (double)sqrtf(x * x + y * y + z * z); }
	__host__ __device__ inline double lenSquared() { return (double)(x * x + y * y + z * z); }
	__host__ __device__ Float3 zeroIfAbove(double a) { return Float3(x * (x < a), y * (y < a), z * (z < a)); }
	__host__ __device__ Float3 zeroIfBelow(double a) { return Float3(x * (x > a), y * (y > a), z * (z > a)); }
	__host__ __device__ Float3 elementwiseModulus(double a) {
		while (x > a)
			x -= a;
		while (y > a)
			y -= a;
		while (z > a)
			z -= a;
		return *this;
	}
	
	__host__ __device__ inline static double getAngle(Float3 v1, Float3 v2) {
		return acos((v1.dot(v2)) / (v1.len() * v2.len()));
	}
	__host__ __device__ static double getAngle(Float3 a, Float3 middle, Float3 b) {
		return getAngle(a - middle, b - middle);
	}


	__host__ __device__ Float3 cross(Float3 a) const { return Float3(y * a.z - z * a.y, z * a.x - x * a.z, x * a.y - y * a.x); }
	__host__ __device__ double dot(Float3 a) const { return (x * a.x + y * a.y + z * a.z); }
	__host__ __device__ Float3 abs() const { return Float3(
		std::abs(x),
		std::abs(y), 
		std::abs(z)
		); }


	__host__ __device__ void print(char c='_') { printf("%c %f %f %f\n", c, x, y, z); }

	__host__ __device__ Float3 rotateAroundOrigin(Float3 pitch_yaw_roll) {	//pitch around x, yaw around z, tilt around y
		// pitch and yaw is relative to global coordinates. Tilt is relative to body direction

		Float3 v = rodriguesRotatation(*this, Float3(1, 0, 0), pitch_yaw_roll.x);

		// Yaw around z
		v = rodriguesRotatation(v, Float3(0, 0, 1), pitch_yaw_roll.y);

		return v;
	}
	__host__ __device__ Float3 rotateAroundVector(Float3 pitch_yaw_roll, Float3 k) {	// k=normal = z-pointing		
		Float3 v = rodriguesRotatation(*this, Float3(1,0,0), pitch_yaw_roll.x);

		v = rodriguesRotatation(v, Float3(0,0,1) , pitch_yaw_roll.y);		
		v = rodriguesRotatation(v, k, pitch_yaw_roll.z);
		return v;
	}
	__host__ __device__ static Float3 rodriguesRotatation(Float3 v, Float3 k, double theta) {
		return v * cos(theta) + k.cross(v) * sin(theta) + k * (k.dot(v)) * (1 - cos(theta));
	}

	__host__ __device__ double at(int index) {
		switch (index) {
		case 0:
			return x;
		case 1:
			return y;
		case 2:
			return z;
		default:
			return -404;
		}
	}

	__host__ __device__ double* placeAt(int index) {
		switch (index) {
		case 0:
			return &x;
		case 1:
			return &y;
		case 2:
			return &z;
		}
	}

	// Not used right now!
	__host__ __device__ static Float3 centerOfMass(Float3* arr_ptr, uint32_t arr_size) {	// Only run before sim, so we can cast to double without slowing sim
		Float3 sum = Float3(0,0,0);
		for (uint32_t i = 0; i < arr_size; i++) {
			sum = sum + arr_ptr[i];
		}
		return sum * (1.f/ arr_size);
	}



	double x = 0, y = 0, z = 0;
};

/*
struct BlockMutex {
	__device__ BlockMutex(){}
	int mutex = 0;

	

	__device__ void lock() {
		while (atomicCAS(mutex, 0, 1) != 0) {}
	}
	__device__ void unlock() {
		atomicExch(mutex, 0);
	}
};
*/

#include <sstream>

using namespace std;










template<typename T>
T* genericMoveToDevice(T* data_ptr, int n_elements) {
	T* gpu_ptr;
	int bytesize = n_elements * sizeof(T);

	cudaMallocManaged(&gpu_ptr, bytesize);
	cudaMemcpy(gpu_ptr, data_ptr, bytesize, cudaMemcpyHostToDevice);
	delete[] data_ptr;

	data_ptr = gpu_ptr;

	printf("Moved %d bytes to device\n", bytesize);
	return gpu_ptr;
}




struct Trajectory {
	Trajectory() {}
	Trajectory(std::string path) {
		positions = new Float3[max_size];
		double buffer[3];

		printf("Path: ");
		cout << path << std::endl;

		string file_contents = readFileIntoString(path);
		istringstream sstream(file_contents);
		string record;

		int counter = 0;

		int step_cnt = 0;
		while (std::getline(sstream, record)) {
			istringstream line(record);

			int dim = 0;
			while (getline(line, record, ';')) {
				buffer[dim++] = stod(record);

				if (dim == 3) {
					positions[counter++] = Float3(buffer);
					dim = 0;
					//positions[counter - 1].print();
				}
			}
			if (step_cnt == 0) {
				n_particles = counter;
			}
			step_cnt++;
		}

		particle_type = new int[n_particles];
		n_steps = step_cnt;
		printf("Loaded trajectory with %d particles and %d steps\n", n_particles, n_steps);
	}
	Trajectory(int size) {
		positions = new Float3[size];
	}

	void moveToDevice() {
		positions = genericMoveToDevice(positions, max_size);
		particle_type = genericMoveToDevice(particle_type, n_particles);
		cudaDeviceSynchronize();
		printf("Success\n");
	}

	string readFileIntoString(const string& path) {
		auto ss = ostringstream{};
		ifstream input_file(path);
		if (!input_file.is_open()) {
			cerr << "Could not open the file - '"
				<< path << "'" << endl;
			exit(EXIT_FAILURE);
		}
		ss << input_file.rdbuf();
		return ss.str();
	}

	int max_size = 100 * 10000;

	Float3* positions;
	int* particle_type;	//0 solvent, 1 compound, 2 virtual

	int n_particles;
	int n_steps;
};

