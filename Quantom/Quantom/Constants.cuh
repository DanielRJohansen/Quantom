#pragma once




#define LIMA_DEBUGMODE

//#define LIMA_SAFERUN		// Use this for?
#define LIMA_VERBOSE









// SImulation specifics

constexpr float BOX_LEN = 10.f;		// Must be > twice the len of largest compound
constexpr float BOX_LEN_HALF = BOX_LEN / 2.f;

const int N_SOLVATE_MOLECULES = 1000;	
const int MAX_RESIDUES_TO_LOAD = 200;








const int MAX_SOLVENTS = 0xFFFF;
constexpr float CUTOFF = 1.f;	//nm/
//const int MAX_ATOM_TYPES = 16;




const int THREADS_PER_SOLVENTBLOCK = 128;	// Must be >= N_SOLVATE_MOLECULES


const int THREADS_PER_COMPOUNDBLOCK = 128; // Must be >= max comp particles



const int SIMULATION_STEPS = 100000;






const int MAX_COMPOUND_PARTICLES = 128;



const int STEPS_PER_NLIST_UPDATE = 200;
const int STEPS_PER_LOGTRANSFER = 20;
//const int STEPS_PER_TRAJTRANSFER = 100;
const int STEPS_PER_THERMOSTAT = 10;
const int STEPS_PER_TRAINDATATRANSFER = 1000;



const int STEPS_PER_RENDER = 50;

const bool ENABLE_BOXTEMP = true;
const bool APPLY_THERMOSTAT = true;										// Switch to using forcefield_host first
const bool PRINT_TEMP = false;


const int MAX_COMPOUNDS = 0xFF;
const int MAX_ATOMS = 1'000'000;
const int MAX_ATOMS_IN_RESIDUE = 32;











// Related to compound bridges
const int COMPOUNDBRIDGES_IN_BUNDLE = 64;
const int MAX_PARTICLES_IN_BRIDGE = 32;
const int MAX_SINGLEBONDS_IN_BRIDGE = 32;
const int MAX_ANGLEBONDS_IN_BRIDGE = 32;
const int MAX_DIHEDRALBONDS_IN_BRIDGE = 32;




