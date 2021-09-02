#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "math.h"



constexpr float PI = 3.14159;


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
	__host__ __device__ inline bool operator == (const Double3 a) const { return (a.x == x && a.y == y && a.z == z); }

	__host__ __device__ inline double len() {return (double)sqrtf(x * x + y * y + z * z); }
	__host__ __device__ Double3 cross(Double3 a) const { return Double3(y * a.z - z * a.y, z * a.x - x * a.z, x * a.y - y * a.x); }

	__host__ __device__ double dot(Double3 a) const { return (x * a.x + y * a.y + z * a.z); }

	__host__ __device__ Double3 rotateAroundOrigin(Double3 pitch_yaw_roll) {	//pitch around x, yaw around z, tilt around y
		// Pitch around x
		Double3 v = rodriguesRotatation(*this, Double3(1, 0, 0), pitch_yaw_roll.x);

		// Yaw around y
		v = rodriguesRotatation(v, Double3(0, 1, 0), pitch_yaw_roll.y);

		// Tilt around itself
		v = rodriguesRotatation(v, v, pitch_yaw_roll.z);
		return v;


		/*float sin_pitch = sin(pitch_yaw_roll.x);
		float cos_pitch = cos(pitch_yaw_roll.x);

		float sin_yaw = sin(pitch_yaw_roll.z);
		float cos_yaw = cos(pitch_yaw_roll.z);

		float sin_tilt = sin(pitch_yaw_roll.y);
		float cos_tilt = cos(pitch_yaw_roll.y);

		//Rotate around x
		float y_x = cos_pitch * y + sin_pitch * z;
		float z_x = -sin_pitch * y + cos_pitch * z;*/

		// Rotate 

	}

	__host__ __device__ Double3 rodriguesRotatation(Double3 v, Double3 k, double theta) {
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
		//default:
			//printf("Double3 error!");
			//exit(-2);
		}
	}

	//Utilities
	//void print()

	double x, y, z;
};









/*
struct cudaBody {
	cudaBody(){}
	cudaBody(SimBody* simbody) : simbody(simbody) {

	}


	float pos[3];
	SimBody* simbody;
};

class KdTree {
	class Node;
	void swap(Node parent, Node child);



public:
	KdTree(){}
	KdTree(Double3 block_center) : block_center(block_center) {}

	void addNode(SimBody* body) {
		double pos[3] = { body->pos.x, body->pos.y,	body->pos.z };

		if (root == NULL)
			root = new Node(body, NULL);
		else
			root->addNode(body, pos);
	}

	void balance() {
		if (root != NULL) {
			root->rootBalance(block_center);

		}
	}

private:
	Node* root = NULL;
	Double3 block_center;


	// Helper functions
	static void swap(Node* parent, Node* child) {
		Node* temp = parent;
		parent = child;
		parent->dim = child->dim;

		child = temp;
		child->dim = temp->dim;
	}

	// Structs

	struct Dim {
		Dim(){}
		Dim(uint8_t val) : val(val) {}
		uint8_t val = 0;
		Dim next() {
			if (val == 2)
				return Dim(0);
			return Dim(val + 1);
		}
	};


	class Node {
	public:
		Node() {}
		Node(SimBody* body, Node* parent, double pos[3]) : body(body), parent(parent) {
			for (int i = 0; i < 3; i++)
				this->pos[i] = pos[i];
			dim = parent->dim.next();
		}
		Dim dim;

		void addNode(SimBody* body, double* pos) {
			if (pos[dim.val] > this->pos[dim.val]) {
				if (right == NULL)
					right = new Node(body, this, pos);
				else
					right->addNode(body, pos);
			}
			else {
				if (left == NULL)
					left = new Node(body, this, pos);
				else
					left->addNode(body, pos);
			}
		}

		void rootBalance(Double3 block_center) {
			double left_dist = 99999;
			double right_dist = 99999;
			double root_dist = (body->pos - block_center).len();
			if (left != NULL)
				left_dist = (left->body->pos - block_center).len();
			if (right != NULL)
				right_dist = (right->body->pos - block_center).len();

			if (left_dist < right_dist) {
				if (left_dist < root_dist)
					swap(left, this);
			}
			else {
				if (right_dist < root_dist)
					swap(right, this);
			}
				


			if (left != NULL)
				left->nodeBalance();
			if (right != NULL)
				right->nodeBalance();
		}
		void nodeBalance();


	private:
		double pos[3];
		SimBody* body;
		Node* parent;
		Node* right = NULL;
		Node* left = NULL;

	};
};

*/