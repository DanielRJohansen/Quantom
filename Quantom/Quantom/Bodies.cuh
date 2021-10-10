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
	//RenderBody* renderbody;

	uint32_t id = UINT32_MAX;

	Float3 pos;	//CoM - nm
	Float3 vel_prev;


	float mass = 0;		// g/mol

	// Cosmetic variables
	float radius = 0;
	uint8_t color[3] = { 255, 100, 0 };

	//unsigned char molecule_type = UNUSED_BODY;			// 255 is unutilized, 0 is water
	bool active = false;
	uint32_t compoundID = UINT32_MAX;
	//uint32_t bondpair_ids[4] = {UINT32_MAX, UINT32_MAX , UINT32_MAX , UINT32_MAX };		// only to avoid intermol forces between bonded atoms
	// these are handled later when we create the compound

	//Float3 charge_unit_vector;
	//float charge_magnitude = 0.f;

	Float3 force;	// J/mol
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

struct CompactParticle {
	CompactParticle(){}
	CompactParticle(Float3 pos, float mass) : pos(pos), mass(mass) {}
	Float3 pos;
	Float3 force;	// Sum all intramol forces here before pushing to box!
	float mass;
};

struct BondPair {	// IDS and indexes are used interchangeably here!
	BondPair(){}
	BondPair(float ref_dist, uint32_t particleindex_a, uint32_t particleindex_b) : 
		reference_dist(ref_dist) {
		atom_indexes[0] = particleindex_a;
		atom_indexes[1] = particleindex_b;
	}

	//uint32_t bond_index;	
	float reference_dist;
	uint32_t atom_indexes[2];	// Relative to the compund - NOT ABSOLUTE INDEX. Used in global table with compunds start-index
};

struct AngleBond {
	float reference_theta;
	uint32_t atom_indexes[3]; // i,j,k angle between i and k
};






// A shell script will automate writing these compounds
const int H2O_PARTICLES = 3;
const int H2O_PAIRBONDS = 2;
const float OH_refdist = 0.095;
struct Compound_H2O {	// Entire molecule for small < 500 atoms molcules, or part of large molecule
	__host__ __device__ Compound_H2O() {};	// {O, H, H}

	__host__ void init(uint32_t startindex_particle, uint32_t compoundID) {
		this->startindex_particle = startindex_particle;
		bondpairs[0] = BondPair(OH_refdist, 0, 1);
		bondpairs[1] = BondPair(OH_refdist, 0, 2);
	}

	/*Compound_H2O operator = (const Compound_H2O a) {
		return Compound_H2O(a->global_particle_table, a.startindex_particle, a.s); 
	}*/
	


	__device__ void loadParticles(Particle* global_particle_table) {

		for (int i = 0; i < n_particles ; i++) {
			particles[i].pos = global_particle_table[startindex_particle + i].pos;
			//particles[i].pos.print();
			//global_particle_table[startindex_particle + i].pos.print();
		}
	}

	uint16_t n_particles = H2O_PARTICLES;
	CompactParticle particles[H2O_PARTICLES];

	uint32_t startindex_particle = 0;
	
	uint16_t n_bondpairs = H2O_PAIRBONDS;
	BondPair bondpairs[H2O_PAIRBONDS];
	
	//Particle* global_particle_table;	// Host address, only usable when creating compound

};

