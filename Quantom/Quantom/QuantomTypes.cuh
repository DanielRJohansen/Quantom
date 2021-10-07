#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "math.h"
#include <iostream>


constexpr float PI = 3.14159;


struct Int3 {
	__host__ __device__ Int3() {}
	__host__ __device__ Int3(int x, int y, int z) : x(x), y(y), z(z) {}


	__host__ __device__ inline Int3 operator + (const Int3 a) const { return Int3(x + a.x, y + a.y, z + a.z); }
	__host__ __device__ inline Int3 operator - (const Int3 a) const { return Int3(x - a.x, y - a.y, z - a.z); }


	int x=0, y=0, z = 0;
};

struct Float3 {
	__host__ __device__ Float3() {}
	__host__ __device__ Float3(float x, float y, float z) : x(x), y(y), z(z) {}

	__host__ __device__ inline Float3 operator * (const float a) const { return Float3(x * a, y * a, z * a); }
	__host__ __device__ inline Float3 operator * (const Float3 a) const { return Float3(x * a.x, y * a.y, z * a.z); }
	__host__ __device__ inline Float3 operator + (const Float3 a) const { return Float3(x + a.x, y + a.y, z + a.z); }
	__host__ __device__ inline Float3 operator - (const Float3 a) const { return Float3(x - a.x, y - a.y, z - a.z); }
	__host__ __device__ inline bool operator == (const Float3 a) const { return (a.x == x && a.y == y && a.z == z); }

	//__host__ __device__ static Float3 norm(Float3 a) { return a * (1.f / a.len()); }	// Remove this at some point..
	__host__ __device__ Float3 norm() {
		float l = len();
		if (l)
			return *this * (1.f / l); }
	__host__ __device__ Float3 square() {return Float3(x * x, y * y, z * z);}
	__host__ __device__ inline float len() {return (float)sqrtf(x * x + y * y + z * z); }
	__host__ __device__ inline float lenSquared() { return (float)(x * x + y * y + z * z); }
	__host__ __device__ Float3 zeroIfAbove(float a) { return Float3(x * (x < a), y * (y < a), z * (z < a)); }
	__host__ __device__ Float3 zeroIfBelow(float a) { return Float3(x * (x > a), y * (y > a), z * (z > a)); }
	__host__ __device__ Float3 elementwiseModulus(float a) {
		while (x > a)
			x -= a;
		while (y > a)
			y -= a;
		while (z > a)
			z -= a;
		return *this;
	}


	__host__ __device__ Float3 cross(Float3 a) const { return Float3(y * a.z - z * a.y, z * a.x - x * a.z, x * a.y - y * a.x); }
	__host__ __device__ float dot(Float3 a) const { return (x * a.x + y * a.y + z * a.z); }
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
	__host__ __device__ static Float3 rodriguesRotatation(Float3 v, Float3 k, float theta) {
		return v * cos(theta) + k.cross(v) * sin(theta) + k * (k.dot(v)) * (1 - cos(theta));
	}

	__host__ __device__ float at(int index) {
		switch (index) {
		case 0:
			return x;
		case 1:
			return y;
		case 2:
			return z;
		default:
			return -1;
		}
	}



	float x = 0, y = 0, z = 0;
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
	KdTree(Float3 block_center) : block_center(block_center) {}

	void addNode(SimBody* body) {
		float pos[3] = { body->pos.x, body->pos.y,	body->pos.z };

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
	Float3 block_center;


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
		Node(SimBody* body, Node* parent, float pos[3]) : body(body), parent(parent) {
			for (int i = 0; i < 3; i++)
				this->pos[i] = pos[i];
			dim = parent->dim.next();
		}
		Dim dim;

		void addNode(SimBody* body, float* pos) {
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

		void rootBalance(Float3 block_center) {
			float left_dist = 99999;
			float right_dist = 99999;
			float root_dist = (body->pos - block_center).len();
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
		float pos[3];
		SimBody* body;
		Node* parent;
		Node* right = NULL;
		Node* left = NULL;

	};
};

*/