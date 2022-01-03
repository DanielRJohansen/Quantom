#pragma once

#include "QuantomTypes.cuh"
#include <SFML/Graphics.hpp>
#include <iostream>

const int MAX_COMPOUND_PARTICLES = 64;














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
	//Float3 force_prev;				// For velocity verlet stormer integration
	Float3 pos_tsub1;				// Must be initiated!
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
const int NEIGHBORLIST_MAX_COMPOUNDS = 32;
const int NEIGHBORLIST_MAX_SOLVENTS = 256;

class NeighborList {
public:
	enum NEIGHBOR_TYPE {COMPOUND, SOLVENT};
	__host__ void init() {
		for (int i = 0; i < NEIGHBORLIST_MAX_COMPOUNDS; i++) {
			neighborcompound_indexes[i] = 0xFFFF;
		}
		for (int i = 0; i < NEIGHBORLIST_MAX_SOLVENTS; i++) {
			neighborsolvent_indexes[i] = 0xFFFF;
		}
	}

	__host__ bool addIndex(int new_index, NEIGHBOR_TYPE nt) {
		switch (nt)
		{
		case NeighborList::COMPOUND:
			if (n_compound_neighbors < NEIGHBORLIST_MAX_COMPOUNDS) {
				neighborcompound_indexes[n_compound_neighbors++] = new_index;
				return true;
			}
			break;
		case NeighborList::SOLVENT:
			if (n_solvent_neighbors < NEIGHBORLIST_MAX_SOLVENTS) {
				neighborsolvent_indexes[n_solvent_neighbors++] = new_index;
				return true;
			}
			break;
		default:
			break;
		}
		printf("Failed!\n");
		return false;
	}

	__device__ void loadMeta(NeighborList* nl_ptr) {	// Called from thread 0
		n_compound_neighbors = nl_ptr->n_compound_neighbors;
		n_solvent_neighbors = nl_ptr->n_solvent_neighbors;
	}
	__device__ void loadData(NeighborList* nl_ptr) {
		if (threadIdx.x < n_compound_neighbors)
			neighborcompound_indexes[threadIdx.x] = nl_ptr->neighborcompound_indexes[threadIdx.x];
		for (int i = threadIdx.x;  i < n_solvent_neighbors; i += blockDim.x) {// Same as THREADS_PER_COMPOUNDBLOCK
			neighborsolvent_indexes[i] = nl_ptr->neighborsolvent_indexes[i];
			i += blockDim.x;	
		}
	}

	uint16_t neighborcompound_indexes[NEIGHBORLIST_MAX_COMPOUNDS];
	uint8_t n_compound_neighbors = 0;
	uint16_t neighborsolvent_indexes[NEIGHBORLIST_MAX_SOLVENTS];
	uint8_t n_solvent_neighbors = 0;
};







// A shell script will automate writing these compounds
const int H2O_PARTICLES = 3;
const int H2O_PAIRBONDS = 2;
const int H2O_ANGLEBONDS = 1;
const double OH_refdist = 0.095;			// nm
const double HOH_refangle = 1.822996;	// radians
const double max_LJ_dist = 1;			// nm

const int MAX_PAIRBONDS = 64;
const int MAX_ANGLEBONDS = 64;
const double CC_refdist = 0.153; // nm
const double CCC_reftheta = 1.953; // nm

struct CompoundState {
	__device__ void setMeta(int n_p) {
		n_particles = n_p;
	}
	__device__ void loadData(CompoundState* state) {
		if (threadIdx.x < n_particles)
			positions[threadIdx.x] = state->positions[threadIdx.x];
	}


	Float3 positions[MAX_COMPOUND_PARTICLES];
	uint8_t n_particles = 0;
};


struct Compound {
	__host__ Compound() {}	// {}


	//---------------------------------------------------------------------------------//
	__host__ Compound(uint32_t index, CompoundState* states_host) {		// Default for testing only!
		this->index = index;

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

	__host__ Float3 getCOM() {
		Float3 com;
		for (int i = 0; i < n_particles; i++)
			com += (particles[i].pos_tsub1 * (1.f / (double) n_particles));
		return com;
	}
	__host__ bool intersects(Compound a) {
		return (a.center_of_mass - center_of_mass).len() < (a.radius + radius + max_LJ_dist);
	}
	//---------------------------------------------------------------------------------//

	__device__ void loadMeta(Compound* compound) {
		n_particles = compound->n_particles;
		n_pairbonds = compound->n_pairbonds;
		n_anglebonds = compound->n_anglebonds;
	}
	__device__ void loadData(Compound* compound) {
		if (threadIdx.x < n_particles) {
			particles[threadIdx.x] = compound->particles[threadIdx.x];
		}
		for (int i = 0; (i * blockDim.x) < n_pairbonds; i++) {
			int index = i * blockDim.x + threadIdx.x;
			if (index < n_pairbonds)
				pairbonds[index] = compound->pairbonds[index];
		}
		for (int i = 0; (i * blockDim.x) < n_anglebonds; i++) {
			int index = i * blockDim.x + threadIdx.x;
			if (index < n_anglebonds)
				anglebonds[index] = compound->anglebonds[index];
		}
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