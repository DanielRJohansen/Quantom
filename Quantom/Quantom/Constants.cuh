#pragma once


#define LIMA_DEBUGMODE

//#define LIMA_SAFERUN		// Use this for?
#define LIMA_VERBOSE


// SImulation specifics
const int N_SOLVATE_MOLECULES = 600;	



const int MAX_COMPOUND_PARTICLES = 128;



const int STEPS_PER_NLIST_UPDATE = 200;
const int STEPS_PER_LOGTRANSFER = 200;
//const int STEPS_PER_TRAJTRANSFER = 100;
const int STEPS_PER_THERMOSTAT = 200;
const int STEPS_PER_TRAINDATATRANSFER = 2000;



const int STEPS_PER_RENDER = 50;


const bool APPLY_THERMOSTAT = true;										// Switch to using forcefield_host first



const int MAX_COMPOUNDS = 0xFF;
const int MAX_ATOMS = 1'000'000;
const int MAX_ATOMS_IN_RESIDUE = 32;











// Related to compound bridges
const int COMPOUNDBRIDGES_IN_BUNDLE = 64;
const int MAX_PARTICLES_IN_BRIDGE = 32;
const int MAX_SINGLEBONDS_IN_BRIDGE = 32;
const int MAX_ANGLEBONDS_IN_BRIDGE = 32;
const int MAX_DIHEDRALBONDS_IN_BRIDGE = 32;




