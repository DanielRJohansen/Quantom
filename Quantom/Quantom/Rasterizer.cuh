#pragma once

#include "QuantomTypes.cuh"
#include "cuda_runtime.h"
#include "Simulation.cuh"


const int RAS_THREADS_PER_BLOCK = 64;


struct RenderAtom {
	Float3 pos;
	//ATOM_TYPE atom_type;
	float mass;
	float radius;
	Int3 color;
	ATOM_TYPE atom_type = ATOM_TYPE::NONE;
};



class Rasterizer {
public:
	Rasterizer() {};
	
	RenderBall* render(Simulation* simulation);

	int actual_n_particles;

private:
	RenderAtom* getAllAtoms(Simulation* simulation);
	void sortAtoms(RenderAtom* atoms, int dim);
	RenderBall* processAtoms(RenderAtom* atoms);


	int n_threadblocks;
	//int actual_n_particles;
	int solvent_offset;
};
