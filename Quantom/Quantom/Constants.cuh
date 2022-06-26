#pragma once




#define LIMA_DEBUGMODE

//#define LIMA_SAFERUN		// Use this for?
//#define LIMA_VERBOSE

//#define ENABLE_WATER


#define ENABLE_DIHEDRALS





// Debugging
const bool print_compound_positions = false;



// SImulation specifics

constexpr float BOX_LEN = 7.f;		// Must be > twice the len of largest compound
constexpr float BOX_LEN_HALF = BOX_LEN / 2.f;

const int N_SOLVATE_MOLECULES = 2000;	


const int NEIGHBORLIST_MAX_COMPOUNDS = 64;
const int NEIGHBORLIST_MAX_SOLVENTS = 1024;
const int SOLVENTS_BATCHSIZE = 64;







const int MAX_SOLVENTS = 0xFFFF;
constexpr float CUTOFF = 0.7f;	//nm/
//const int MAX_ATOM_TYPES = 16;



// OPTIMIZATION PARAMETERS
const int MAX_COMPOUND_PARTICLES = 64;



const int THREADS_PER_SOLVENTBLOCK = 128;	
const int THREADS_PER_COMPOUNDBLOCK = MAX_COMPOUND_PARTICLES; // Must be >= max comp particles


const int SIMULATION_STEPS = 1;




const int STEPS_PER_NLIST_UPDATE = 200;



const int STEPS_PER_RENDER = 50;


// THERMOSTAT PARAMETERS
const int STEPS_PER_LOGTRANSFER = 40;
//const int STEPS_PER_TRAJTRANSFER = 100;
const int STEPS_PER_THERMOSTAT = 40;
const bool ENABLE_BOXTEMP = true;
const bool APPLY_THERMOSTAT = true;										// Switch to using forcefield_host first
const bool PRINT_TEMP = false;


const int MAX_COMPOUNDS = 0xFF;
const int MAX_ATOMS = 1'000'000;
const int MAX_ATOMS_IN_RESIDUE = 32;






// Logging constants
const int N_DATAGAN_VALUES = 3;
const int STEPS_PER_TRAINDATATRANSFER = 100;








// Related to compound bridges
const int COMPOUNDBRIDGES_IN_BUNDLE = 64;
const int MAX_PARTICLES_IN_BRIDGE = 32;
const int MAX_SINGLEBONDS_IN_BRIDGE = 32;
const int MAX_ANGLEBONDS_IN_BRIDGE = 32;
const int MAX_DIHEDRALBONDS_IN_BRIDGE = 32;




