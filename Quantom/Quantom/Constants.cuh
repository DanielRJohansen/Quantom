#pragma once




#define LIMA_DEBUGMODE

//#define LIMA_SAFERUN		// Use this for?
#define LIMA_VERBOSE

#define ENABLE_WATER







// SImulation specifics

constexpr float BOX_LEN = 8.f;		// Must be > twice the len of largest compound
constexpr float BOX_LEN_HALF = BOX_LEN / 2.f;

const int N_SOLVATE_MOLECULES = 2000;	








const int MAX_SOLVENTS = 0xFFFF;
constexpr float CUTOFF = 1.f;	//nm/
//const int MAX_ATOM_TYPES = 16;



// OPTIMIZATION PARAMETERS
const int THREADS_PER_SOLVENTBLOCK = 128;	// Must be >= N_SOLVATE_MOLECULES
const int THREADS_PER_COMPOUNDBLOCK = 128; // Must be >= max comp particles
const int MAX_COMPOUND_PARTICLES = 128;


const int SIMULATION_STEPS = 10000;




const int STEPS_PER_NLIST_UPDATE = 200;



const int STEPS_PER_RENDER = 100;


// THERMOSTAT PARAMETERS
const int STEPS_PER_LOGTRANSFER = 10;
//const int STEPS_PER_TRAJTRANSFER = 100;
const int STEPS_PER_THERMOSTAT = 10;
const bool ENABLE_BOXTEMP = true;
const bool APPLY_THERMOSTAT = true;										// Switch to using forcefield_host first
const bool PRINT_TEMP = false;


const int MAX_COMPOUNDS = 0xFF;
const int MAX_ATOMS = 1'000'000;
const int MAX_ATOMS_IN_RESIDUE = 32;






// Logging constants
const int N_DATAGAN_VALUES = 6;
const int STEPS_PER_TRAINDATATRANSFER = 100;








// Related to compound bridges
const int COMPOUNDBRIDGES_IN_BUNDLE = 64;
const int MAX_PARTICLES_IN_BRIDGE = 32;
const int MAX_SINGLEBONDS_IN_BRIDGE = 32;
const int MAX_ANGLEBONDS_IN_BRIDGE = 32;
const int MAX_DIHEDRALBONDS_IN_BRIDGE = 32;




