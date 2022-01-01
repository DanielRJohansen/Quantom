#pragma once

#include "QuantomTypes.cuh"
#include <SFML/Graphics.hpp>
#include <iostream>

//enum Atom{Oxygen};

struct Atom {
	__host__ __device__ Atom() {}
	__host__ Atom(Float3 pos, double r, double mass, uint8_t c[3]) : pos(pos), radius(r), mass(mass) {
		mass *= 0.001f;	// Convert to kg/mol
		for (int i = 0; i < 3; i++) {
			color[i] = c[i];
		}
	}
	Float3 pos;	// Relative	to CoM, and (0,0,0) rotation
	double radius;	// in fm?
	double mass;		// in kg/mol
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
	double radii[3];
};

constexpr double BODY_RADIUS = 0.2;		// CRITICAL VALUE!
constexpr unsigned char UNUSED_BODY = 255;
struct CompactParticle {	// Contains information only needed by the Ownerkernel
	CompactParticle() {}	
	CompactParticle(double mass, Float3 pos_sub1) : mass(mass), pos_tsub1(pos_sub1)  {}
	//Float3 vel;						// nm/ns
	//Float3 acc;						// nm/ns^2
	//Float3 force_prev;				// For velocity verlet stormer integration
	Float3 pos_tsub1;				// Must be initiated!
	//double pot_E_prev = 999999999;			// MUST BE INITIATED BY CALCULATION, OR IT WILL FUCK SHIT UP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	double mass;								// g/mol
};







struct Solvent {
	Solvent() {}
	Solvent(Float3 pos, Float3 pos_tsub1) : pos(pos), pos_tsub1(pos_tsub1) {}

	Float3 pos;
	Float3 pos_tsub1;
};






//--------------------------- THE FOLLOWING IS FOR HANDLING INTRAMOLECULAR FORCES ---------------------------//

struct PairBond {	// IDS and indexes are used interchangeably here!
	PairBond(){}
	PairBond(double ref_dist, uint32_t particleindex_a, uint32_t particleindex_b) : 
		reference_dist(ref_dist) {
		atom_indexes[0] = particleindex_a;
		atom_indexes[1] = particleindex_b;
	}

	//uint32_t bond_index;	
	double reference_dist;
	uint32_t atom_indexes[2];	// Relative to the compund - NOT ABSOLUTE INDEX. Used in global table with compunds start-index
};

struct AngleBond {
	AngleBond() {}
	AngleBond(double ref_t, uint32_t particleindex_l, uint32_t particleindex_m, uint32_t particleindex_r) :
		reference_angle(ref_t) {
		atom_indexes[0] = particleindex_l;
		atom_indexes[1] = particleindex_m;
		atom_indexes[2] = particleindex_r;
	}

	double reference_angle;
	uint32_t atom_indexes[3]; // i,j,k angle between i and k
};




// ------------------------------------------------- COMPOUNDS ------------------------------------------------- //

class NeighborList {
public:
	__host__ void init(uint16_t* index_buffer, int max_indexes) {
		for (int i = 0; i < max_indexes; i++) {
			index_buffer[i] = 0xFFFF;
		}
	}
	__host__ void addIndex(int new_index, uint16_t* index_buffer, int max_indexes, uint8_t* neighbor_cnt) {
		for (int i = 0; i < max_indexes; i++) {
			if (index_buffer[i] == 0xFFFF) {
				index_buffer[i] = new_index;
				(*neighbor_cnt)++;
				return;
			}
		}
	}
};
class CompoundNeighborList : public NeighborList {	// Both compounds and solvents have one of these
public:
	CompoundNeighborList() {
		init(neighborcompound_indexes, 32);
	}
	__host__ void addIndex(int new_index) {
		NeighborList::addIndex(new_index, neighborcompound_indexes, 32, &n_neighbors);
	}
	uint16_t neighborcompound_indexes[32]; // For now so we dont need to track, adjust to 16 or 32 later!!
	uint8_t n_neighbors = 0;					// adjust too?
};

class SolventNeighborList : public NeighborList{						// Both compounds and solvents have one of these
public:
	SolventNeighborList() {
		init(neighborsolvent_indexes, 256);
	}
	__host__ void addIndex(int new_index) {
		NeighborList::addIndex(new_index, neighborsolvent_indexes, 256, &n_neighbors);
	}
	uint16_t neighborsolvent_indexes[256]; 
	uint8_t n_neighbors = 0;
};








// A shell script will automate writing these compounds
const int H2O_PARTICLES = 3;
const int H2O_PAIRBONDS = 2;
const int H2O_ANGLEBONDS = 1;
const double OH_refdist = 0.095;			// nm
const double HOH_refangle = 1.822996;	// radians
const double max_LJ_dist = 1;			// nm

const int MAX_COMPOUND_PARTICLES = 64;
const int MAX_PAIRBONDS = 32;
const int MAX_ANGLEBONDS = 32;
const double CC_refdist = 0.153; // nm
const double CCC_reftheta = 1.953; // nm

struct CompoundState {
	Float3 positions[MAX_COMPOUND_PARTICLES];
	uint8_t n_particles = 0;
};

struct Compound {
	__host__ Compound() {}	// {}


	//---------------------------------------------------------------------------------//
	__host__ Compound(uint32_t index, CompoundState* states_host) {
		this->index = index;
		//compound_neighborlist_ptr = neighborlist_device;
		//compound_state_ptr = state_device;

		pairbonds[0] = PairBond(CC_refdist, 0, 1);
		n_pairbonds++;

		pairbonds[1] = PairBond(CC_refdist, 0, 2);
		n_pairbonds++;

		anglebonds[0] = AngleBond(CCC_reftheta, 1, 0, 2);
		n_anglebonds++;

		for (uint32_t i = 0; i < n_particles; i++)
			center_of_mass = center_of_mass + states_host->positions[i];
		center_of_mass = center_of_mass * (1.f / states_host->n_particles);

		radius = pairbonds[0].reference_dist * states_host->n_particles;					// TODO: Shitty estimate, do better later
	};
	
	__host__ void init(uint32_t index) {	// Only call this if the compound has already been assigned particles & bonds
		this->index = index;
		center_of_mass = getCOM();
		//printf("")
		radius = pairbonds[0].reference_dist * n_particles * 0.5f;
		center_of_mass.print('C');
		printf("Radius %f\n", radius);
	}

	//---------------------------------------------------------------------------------//

	__host__ Float3 getCOM() {
		Float3 com;
		for (int i = 0; i < n_particles; i++)
			com += (particles[i].pos_tsub1 * (1.f / (double) n_particles));
		return com;
	}
	__host__ bool intersects(Compound a) {
		return (a.center_of_mass - center_of_mass).len() < (a.radius + radius + max_LJ_dist);
	}

	uint32_t index;										// Is this necessary
	//CompoundState* compound_state_ptr;
	//CompoundNeighborList* compound_neighborlist_ptr;

	uint8_t n_particles = 0;					// MAX 256 particles!!!!0
	CompactParticle particles[MAX_COMPOUND_PARTICLES];

	Float3 center_of_mass = Float3(0, 0, 0);
	double radius = 0;

	uint16_t n_pairbonds = 0;
	PairBond pairbonds[MAX_PAIRBONDS];

	uint16_t n_anglebonds = 0;
	AngleBond anglebonds[MAX_ANGLEBONDS];
};


struct Molecule1 {
	int n_compounds = 0;
	Compound* compounds;
};