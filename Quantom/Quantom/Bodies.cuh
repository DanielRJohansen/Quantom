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

struct RenderBody {

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

constexpr float BODY_RADIUS = 0.2;		// CRITICAL VALUE!
constexpr unsigned char UNUSED_BODY = 255;

struct Particle {
	__host__ __device__ Particle() {}
	__host__ Particle(uint32_t id, Float3 pos, Float3 vel_prev, float mass, uint32_t compoundID) :
		id(id), pos(pos), vel_prev(vel_prev), mass(mass), compoundID(compoundID) {
		active = true;
	}

	uint32_t id = UINT32_MAX;

	Float3 pos;	//CoM - nm
	Float3 vel_prev;
	Float3 force;	// J/mol

	float mass = 0;		// g/mol

	// Cosmetic variables
	float radius = 0;
	uint8_t color[3] = { 255, 100, 0 };
	bool active = false;
	uint32_t compoundID = UINT32_MAX;

};

struct CompactParticle {
	CompactParticle() {}
	CompactParticle(Float3 pos, float mass) : pos(pos) {}
	Float3 pos;
	Float3 vel_prev;
	float pot_E_prev = 0;
	float mass;
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
		reference_theta(ref_t) {
		atom_indexes[0] = particleindex_l;
		atom_indexes[1] = particleindex_m;
		atom_indexes[2] = particleindex_r;
	}

	float reference_theta;
	uint32_t atom_indexes[3]; // i,j,k angle between i and k
};




// ------------------------------------------------- COMPOUNDS ------------------------------------------------- //

struct CompoundNeighborInfo {
	uint32_t neighborcompound_indexes[256]; // For now so we dont need to track, adjust to 16 or 32 later!!
	uint8_t n_neighbors;					// adjust too?
};
struct CompoundState {
	CompactParticle particles[128];
	uint8_t particle_cnt;
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
	__host__ Compound_H2O(uint32_t index, CompoundNeighborInfo* info, CompoundState* state) {
		this->index = index;
		compound_neighborinfo_ptr = info;
		compound_state_ptr = state;

		pairbonds[0] = PairBond(OH_refdist, 0, 1);
		pairbonds[1] = PairBond(OH_refdist, 0, 2);

		anglebonds[0] = AngleBond(HOH_refangle, 1, 0, 2);

		for (uint32_t i = 0; i < n_particles; i++)
			center_of_mass = center_of_mass + particles[i].pos;
		center_of_mass = center_of_mass * (1.f / n_particles);

		radius = pairbonds[0].reference_dist * n_particles;					// TODO: Shitty estimate, do better later
	};	
	
	__host__ bool intersects(Compound_H2O a) {
		return (a.center_of_mass - center_of_mass).len() < (a.radius + radius + max_LJ_dist);
	}

	uint32_t index;
	CompoundState* compound_state_ptr;
	CompoundNeighborInfo* compound_neighborinfo_ptr;

	uint8_t n_particles = H2O_PARTICLES;
	CompactParticle particles[H2O_PARTICLES];

	Float3 center_of_mass = Float3(0,0,0);
	float radius;

	uint16_t n_pairbonds = H2O_PAIRBONDS;
	PairBond pairbonds[H2O_PAIRBONDS];
	
	uint16_t n_anglebonds = H2O_ANGLEBONDS;
	AngleBond anglebonds[H2O_ANGLEBONDS];
};



