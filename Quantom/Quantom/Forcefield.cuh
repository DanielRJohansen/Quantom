#pragma once

#include "Bodies.cuh"
#include "Simulation.cuh"



#define ATOMTYPE_C 0
#define ATOMTYPE_O 1
#define ATOMTYPE_N 2
#define ATOMTYPE_H 3


struct ParticleParameters {
	const float mass = -1;

	//Nonbonded
	const float sigma = -1;
	const 	float epsilon = -1;

	// Bonded
	const float b0 = -1;		// ref_dist
	const float kb = -1;

	// Angles
	const float theta_0 = -1;
	const float ktheta = -1;
	//const float ub_0;		// ???
	//const float kub;		// ???

	// Dihedrals
	const float phi_0 = -1;		// Ref torsion. Is always 0 or 180, could be a binary
	const float k_phi = -1;		// 
};

struct PPMaker {
	PPMaker() {}
	PPMaker(char atom, float m, float s, float e) : mass(m), sigma(s), epsilon(e) {}
	PPMaker(char atom, float m, float s, float e, float b0, float kb, float t0, float kt, float p0, float kp) : mass(m), sigma(s), epsilon(e),
		b0(b0), kb(kb), theta_0(t0), ktheta(kt), phi_0(p0), k_phi(kp) {}
	//uint8_t id;
	const float mass = -1;

	//Nonbonded
	const float sigma = -1;
	const float epsilon = -1;

	// Bonded
	const float b0 = -1;		// ref_dist
	const float kb = -1;

	// Angles
	const float theta_0 = -1;
	const float ktheta = -1;
	//const float ub_0;		// ???
	//const float kub;		// ???

	// Dihedrals
	const float phi_0 = -1;		// Ref torsion. Is always 0 or 180, could be a binary
	const float k_phi = -1;		// 
};




struct ForceField {
	const ParticleParameters particle_parameters[MAX_ATOM_TYPES];
};

class ForceFieldMaker {
public:
	ForceFieldMaker() {
		make('C' , 0, 12.011f, 0.356359487256f, 0.46024f);
		/*forcefield.particle_parameters[0] = ParticleParameters('C', 12.011f, 0.356359487256f, 0.46024f);
		forcefield.particle_parameters[1] = ParticleParameters('O', 12.011f, 0.356359487256f, 0.46024f);
		forcefield.particle_parameters[2] = ParticleParameters('N', 12.011f, 0.356359487256f, 0.46024f);
		forcefield.particle_parameters[3] = ParticleParameters('H', 12.011f, 0.356359487256f, 0.46024f);
		*/
	}

	void make(char atom, int index, float mass, float s, float e) {
		//forcefield.particle_parameters[index].mass = mass;
		//forcefield.particle_parameters[index].sigma = s;
		//forcefield.particle_parameters[index].epsilon = e;

	}

	ForceField getForcefield() {
		return forcefield;
	}

	int atomTypeToIndex(char atom) {
		if (atom == 'C')
			return 0;
		if (atom == 'O')
			return 1;
		if (atom == 'N')
			return 2;
		if (atom == 'H')
			return 3;
		printf("Unable to find atom %c\n", atom);
		exit(1);
	}

private:
	ForceField forcefield;
};

