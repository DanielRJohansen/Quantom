#pragma once

#include "Bodies.cuh"
//#include "Simulation.cuh"
#include "Filehandling.h"
#include <string>
#include <vector>

using std::string;


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


	void buildForcefield();
	PairBond* getBondType(int id1, int id2);
	AngleBond* getAngleType(int id1, int id2, int id3);



	void make(char atom, int index, float m, float s, float e) {
		forcefield.particle_parameters[index].mass = m;
		forcefield.particle_parameters[index].sigma = s;
		forcefield.particle_parameters[index].epsilon = e;
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
	ForceField forcefield1;


	NBAtomtype* nb_atomtypes;
	int n_nb_atomtypes = 0;

	PairBond* topol_bonds;
	int n_topol_bonds = 0;
	AngleBond* topol_angles;
	int n_topol_angles = 0;
	DihedralBond* topol_dihedrals;
	int n_topol_dihedrals = 0; 



	//enum STATE { INACTIVE, NB_ATOMTYPE_MAPPINGS, FF_NONBONDED, NB_ATOMTYPES, PAIRTYPES };
	enum STATE { INACTIVE, FF_NONBONDED, NB_ATOMTYPES, BONDS, ANGLES, DIHEDRALS };
	static STATE setState(string s, STATE current_state) {
		if (s == "ff_nonbonded")
			return FF_NONBONDED;
		if (s == "atoms")
			return NB_ATOMTYPES;
		if (s == "bonds")
			return BONDS;
		if (s == "angles")
			return ANGLES;
		if (s == "dihedrals")
			return DIHEDRALS;
		return current_state;
	}

	bool newParseTitle(vector<string> row) {
		return (row[0][0] == '#');
	}

	//StringMap parseNBAtomtypeMaps(vector<vector<string>> forcefield_rows) {}


	
	NBAtomtype* parseAtomTypes(vector<vector<string>> summary_rows);

	int* parseAtomTypeIDs(vector<vector<string>> forcefield_rows);

	PairBond* parseBonds(vector<vector<string>> forcefield_rows);
	AngleBond* parseAngles(vector<vector<string>> forcefield_rows);
	DihedralBond* parseDihedrals(vector<vector<string>> forcefield_rows);
	
	void loadAtomypesIntoForcefield();

};






