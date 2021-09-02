#pragma once

#include "QuantomTypes.cuh"
#include <SFML/Graphics.hpp>
#include <iostream>

//enum Atom{Oxygen};

struct Atom {
	__host__ __device__ Atom() {}
	__host__ Atom(Double3 pos, double r, double mass, uint8_t c[3]) : pos(pos), radius(r), mass(mass) {
		
		for (int i = 0; i < 3; i++) {
			color[i] = c[i];
		}
	}
	Double3 pos;	// Relative	to CoM, and (0,0,0) rotation
	double radius;	// in fm?
	double mass;
	uint8_t color[3] = { 0,100,0 };
};

struct RenderBody {

};

struct Molecule {
	Molecule();
	int n_atoms;
	Atom* atoms;
	Double3 CoM;	// Relative		

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

constexpr double BODY_RADIUS = 0.6;

struct SimBody {
	//RenderBody* renderbody;



	Double3 pos;	//CoM - nm
	Double3 vel;	// nm/sec

	Double3 rotation;	//radians
	Double3 rot_vel;	//radians/sec
	//double mass;
	//double radius = 0.05;
	char molecule_type = 0;

	//Double3 charge_unit_vector;
	double charge_magnitude;

};





struct Compound {	//Or molecule
	Compound() {};
	Compound(int n_bodies) : n_bodies(n_bodies) {
		bodies = new SimBody[n_bodies]; 
	}

	bool non_bonded = true;

	int n_bodies = 0;
	SimBody* bodies;
};


struct MoleculeLibrary {

	MoleculeLibrary(){
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