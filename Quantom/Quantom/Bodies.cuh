#pragma once

#include "QuantomTypes.cuh"
#include <SFML/Graphics.hpp>
#include <iostream>

//enum Atom{Oxygen};

struct Atom {
	__host__ __device__ Atom() {}
	__host__ Atom(Float3 pos, float r, float mass, uint8_t c[3]) : pos(pos), radius(r), mass(mass) {
		
		for (int i = 0; i < 3; i++) {
			color[i] = c[i];
		}
	}
	Float3 pos;	// Relative	to CoM, and (0,0,0) rotation
	float radius;	// in fm?
	float mass;
	uint8_t color[3] = { 0,100,0 };
};


struct Molecule {
	Molecule();
	int n_atoms;
	Atom* atoms;
	Float3 CoM;	// Relative		

	void moveToDevice() {
		Atom* atoms_temp;
		int bytesize = n_atoms * sizeof(Atom);
		cudaMallocManaged(&atoms_temp, bytesize);
		cudaMemcpy(atoms_temp, atoms, bytesize, cudaMemcpyHostToDevice);
		cudaDeviceSynchronize();
		delete atoms;
		atoms = atoms_temp;
	}
};

struct MoleculeLibrary {

	MoleculeLibrary() {
		molecules = new Molecule[1];
		n_molecules++;


		moveToDevice();
	}

	void moveToDevice() {
		for (int i = 0; i < n_molecules; i++)
			molecules[i].moveToDevice();

		Molecule* temp;
		int bytesize = n_molecules * sizeof(Molecule);
		cudaMallocManaged(&temp, bytesize);
		cudaMemcpy(temp, molecules, bytesize, cudaMemcpyHostToDevice);
		cudaDeviceSynchronize();
		delete molecules;
		molecules = temp;

	}

	int n_molecules = 0;
	Molecule* molecules;

};












struct RenderMolecule {	// Just temporary, dont know howq to properly implement functionality for rendering.
	uint8_t colors[3][3];
	float radii[3];
};

constexpr float BODY_RADIUS = 0.2;		// CRITICAL VALUE!
constexpr unsigned char UNUSED_BODY = 255;


struct CompactParticle {	// Contains information only needed by the Ownerkernel
	CompactParticle() {}	
	CompactParticle(float mass, Float3 initial_vel) : mass(mass), vel_prev(initial_vel)  {}
	Float3 vel_prev;
	float pot_E_prev = 999999999;			// MUST BE INITIATED BY CALCULATION, OR IT WILL FUCK SHIT UP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	float mass;
};








//--------------------------- THE FOLLOWING IS FOR HANDLING INTRAMOLECULAR FORCES ---------------------------//

struct PairBond {	// IDS and indexes are used interchangeably here!
	PairBond(){}
	PairBond(float ref_dist, uint32_t particleindex_a, uint32_t particleindex_b) : 
		reference_dist(ref_dist) {
		atom_indexes[0] = particleindex_a;
		atom_indexes[1] = particleindex_b;
	}

	//uint32_t bond_index;	
	float reference_dist;
	uint32_t atom_indexes[2];	// Relative to the compund - NOT ABSOLUTE INDEX. Used in global table with compunds start-index
};

struct AngleBond {
	AngleBond() {}
	AngleBond(float ref_t, uint32_t particleindex_l, uint32_t particleindex_m, uint32_t particleindex_r) :
		reference_angle(ref_t) {
		atom_indexes[0] = particleindex_l;
		atom_indexes[1] = particleindex_m;
		atom_indexes[2] = particleindex_r;
	}

	float reference_angle;
	uint32_t atom_indexes[3]; // i,j,k angle between i and k
};




// ------------------------------------------------- COMPOUNDS ------------------------------------------------- //

struct CompoundNeighborList {
	uint32_t neighborcompound_indexes[256]; // For now so we dont need to track, adjust to 16 or 32 later!!
	uint8_t n_neighbors = 0;					// adjust too?
};
struct CompoundState {
	Float3 positions[128];
	uint8_t n_particles = 0;
};



// A shell script will automate writing these compounds
const int H2O_PARTICLES = 3;
const int H2O_PAIRBONDS = 2;
const int H2O_ANGLEBONDS = 1;
const float OH_refdist = 0.095;			// nm
const float HOH_refangle = 1.822996;	// radians
const float max_LJ_dist = 1;			// nm

struct Compound_H2O {			// Entire molecule for small < 500 atoms molcules, or part of large molecule
	__host__ Compound_H2O() {}	// {O, H, H}
	__host__ Compound_H2O(uint32_t index, CompoundNeighborList* neighborlist_device, CompoundState* state_device,
	CompoundState* states_host) {
		this->index = index;
		compound_neighborlist_ptr = neighborlist_device;
		compound_state_ptr = state_device;

		pairbonds[0] = PairBond(OH_refdist, 0, 1);
		pairbonds[1] = PairBond(OH_refdist, 0, 2);

		anglebonds[0] = AngleBond(HOH_refangle, 1, 0, 2);

		for (uint32_t i = 0; i < n_particles; i++)
			center_of_mass = center_of_mass + states_host->positions[i];
		center_of_mass = center_of_mass * (1.f / n_particles);

		radius = pairbonds[0].reference_dist * n_particles;					// TODO: Shitty estimate, do better later
	};	
	
	__host__ bool intersects(Compound_H2O a) {
		return (a.center_of_mass - center_of_mass).len() < (a.radius + radius + max_LJ_dist);
	}

	uint32_t index;										// Is this necessary
	CompoundState* compound_state_ptr;
	CompoundNeighborList* compound_neighborlist_ptr;

	uint8_t n_particles = H2O_PARTICLES;
	CompactParticle particles[H2O_PARTICLES];

	Float3 center_of_mass = Float3(0,0,0);
	float radius;

	uint16_t n_pairbonds = H2O_PAIRBONDS;
	PairBond pairbonds[H2O_PAIRBONDS];
	
	uint16_t n_anglebonds = H2O_ANGLEBONDS;
	AngleBond anglebonds[H2O_ANGLEBONDS];
};



