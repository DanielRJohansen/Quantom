#pragma once



//#define LIMA_SAFERUN		// Use this for?
#define LIMA_VERBOSE


// SImulation specifics
const int N_SOLVATE_MOLECULES = 600;	




const int MAX_COMPOUND_PARTICLES = 64;



const int STEPS_PER_NLIST_UPDATE = 100;
const int STEPS_PER_LOGTRANSFER = 500;
//const int STEPS_PER_TRAJTRANSFER = 100;
const int STEPS_PER_THERMOSTAT = 500;
const int STEPS_PER_TRAINDATATRANSFER = 500;



const int STEPS_PER_RENDER = 500;


const bool APPLY_THERMOSTAT = false;										// Switch to using forcefield_host first



const int MAX_COMPOUNDS = 0xFF;
const int MAX_ATOMS = 1'000'000;
const int MAX_ATOMS_IN_RESIDUE = 16;











// Related to compound bridges
const int COMPOUNDBRIDGES_IN_BUNDLE = 64;
const int MAX_PARTICLES_IN_BRIDGE = 32;
const int MAX_SINGLEBONDS_IN_BRIDGE = 32;
const int MAX_ANGLEBONDS_IN_BRIDGE = 32;
const int MAX_DIHEDRALBONDS_IN_BRIDGE = 32;




