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

