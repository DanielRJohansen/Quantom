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
	__host__ Particle(uint32_t id, Float3 pos, Float3 vel_prev, float mass) :
		id(id), pos(pos), vel_prev(vel_prev), mass(mass)
	{
		active = true;
	}
	//RenderBody* renderbody;

	uint32_t id = 0;

	Float3 pos;	//CoM - nm
	Float3 vel_prev;


	float mass = 0;
	//float radius = 0.05;
	//unsigned char molecule_type = UNUSED_BODY;			// 255 is unutilized, 0 is water
	bool active = false;
	uint32_t bondpair_ids[4];		// only to avoid intermol forces between bonded atoms
	// these are handled later when we create the compound

	//Float3 charge_unit_vector;
	//float charge_magnitude = 0.f;

	Float3 force;
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
	float mass;
};

struct PairBond {	// IDS and indexes are used interchangeably here!
	PairBond(){}
	PairBond(float ref_dist, uint32_t particleindex_a, uint32_t particleindex_b, uint32_t bond_index) : 
		reference_dist(ref_dist), bond_index(bond_index) {
		atom_indexes[0] = particleindex_a;
		atom_indexes[1] = particleindex_b;
	}

	uint32_t bond_index;	
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
struct Compound_H2O {	// Entire molecule for small < 500 atoms molcules, or part of large molecule
	Compound_H2O() {};	// {O, H, H}
	Compound_H2O(Particle* global_particle_table, uint32_t startindex_particle, uint32_t startindex_bond) : 
		global_particle_table(global_particle_table), startindex_particle(startindex_particle) {
		pairbonds[0] = PairBond(0.095, 0, 1, startindex_bond++);
		pairbonds[1] = PairBond(0.095, 0, 2, startindex_bond++);

		uint8_t particle_bond_count[H2O_PARTICLES] = { 0 };
		for (int i = 0; i < H2O_PAIRBONDS; i++) {	
			for (int j = 0; j < 2; j++) {		// Iterate over both particles in bond
				uint32_t rel_p_index = pairbonds[i].atom_indexes[j];
				uint32_t abs_p_index = startindex_particle + rel_p_index;
				uint8_t p_bond_cnt = particle_bond_count[rel_p_index];
				global_particle_table[abs_p_index].bondpair_ids[p_bond_cnt] = pairbonds[i].bond_index;
				particle_bond_count[rel_p_index]++;
			}
		}
	}

	Compound_H2O operator = (const Compound_H2O a) { return Compound_H2O(&a->global_particle_table, a.startindex_particle, a.s); }
	void loadParticles(Particle* global_particle_table) {
		for (int i = 0; i < n_particles ; i++) {
			particles[i].pos = global_particle_table[startindex_particle + i].pos;
		}
	}

	const uint16_t n_particles = H2O_PARTICLES;
	CompactParticle particles[H2O_PARTICLES];

	uint32_t startindex_particle = 0;
	
	const uint16_t n_pairbonds = H2O_PAIRBONDS;
	PairBond pairbonds[H2O_PAIRBONDS];
	
	Particle* global_particle_table;	// Host address, only usable when creating compound

};

