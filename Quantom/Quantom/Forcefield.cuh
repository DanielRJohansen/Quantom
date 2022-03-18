#pragma once

//#include "Bodies.cuh"
//#include "Simulation.cuh"

const int MAX_ATOM_TYPES = 16;


#define ATOMTYPE_SOL 0
#define ATOMTYPE_C 1
#define ATOMTYPE_O 2
#define ATOMTYPE_N 3
#define ATOMTYPE_H 4
#define ATOMTYPE_P 5
#define ATOMTYPE_S 6



struct ParticleParameters {

	float mass = -1;		//[kg/mol]

	//Nonbonded
	float sigma = -1;
	float epsilon = -1;

	// Bonded
	float b0 = -1;		// ref_dist
	float kb = -1;

	// Angles
	float theta_0 = -1;
	float ktheta = -1;
	//const float ub_0;		// ???
	//const float kub;		// ???

	// Dihedrals
	float phi_0 = -1;		// Ref torsion. Is always 0 or 180, could be a binary
	float k_phi = -1;		// 
};


struct ForceField {
	ParticleParameters particle_parameters[MAX_ATOM_TYPES];
};



class ForceFieldMaker {
public:
	ForceFieldMaker() {
		make('W', 0, 15.999000 + 2 * 1.008000, 0.302905564168 + 2 * 0.040001352445, 0.50208f + 2 * 0.19246);
		//make('C', 1, 12.011f, 0.356359487256f, 0.46024f, 0.13350000f, 502080.00f);
		make('C', 1, 12.011f, 0.356359487256f, 0.46024f);
		make('O', 2, 15.999000, 0.302905564168f, 0.50208f);
		make('N', 3, 14.007000, 0.329632525712, 0.83680);
		make('H', 4, 1.008000, 0.040001352445, 0.19246);
		make('P', 5, 30.974000, 0.383086448800, 2.44764);
		make('S', 6, 32.060000, 0.356359487256,  1.88280);
		/*forcefield.particle_parameters[0] = ParticleParameters('C', 12.011f, 0.356359487256f, 0.46024f);
		forcefield.particle_parameters[1] = ParticleParameters('O', 12.011f, 0.356359487256f, 0.46024f);
		forcefield.particle_parameters[2] = ParticleParameters('N', 12.011f, 0.356359487256f, 0.46024f);
		forcefield.particle_parameters[3] = ParticleParameters('H', 12.011f, 0.356359487256f, 0.46024f);
		*/

		for (int i = 0; i < MAX_ATOM_TYPES; i++) {
			forcefield.particle_parameters[i].mass /= 1000.f;		// Convert g/moæ to kg/mol
			forcefield.particle_parameters[i].epsilon *= 1000.f;		// convert kJ/mol to J/mol
		}
			
		printf("Forcefield size: %d bytes\n", sizeof(ForceField));
	}

	void make(char atom, int index, float m, float s, float e) {
		forcefield.particle_parameters[index].mass = m;
		forcefield.particle_parameters[index].sigma = s;
		forcefield.particle_parameters[index].epsilon = e;
	}

	void make(char atom, int index, float m, float s, float e, float b0, float kb) {
		forcefield.particle_parameters[index].mass = m;
		forcefield.particle_parameters[index].sigma = s;
		forcefield.particle_parameters[index].epsilon = e;
		forcefield.particle_parameters[index].b0 = b0;
		forcefield.particle_parameters[index].kb = kb;
	}

	ForceField getForcefield() {
		return forcefield;
	}

	int atomTypeToIndex(char atom) {
		if (atom == 'C')
			return 1;
		if (atom == 'O')
			return 2;
		if (atom == 'N')
			return 3;
		if (atom == 'H')
			return 4;
		if (atom == 'P')
			return 5;
		if (atom == 'S')
			return 6;
		printf("Unable to find atom %c\n", atom);
		exit(1);
	}

private:
	ForceField forcefield;
};

