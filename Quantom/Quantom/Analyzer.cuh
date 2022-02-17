#include <iostream>
#include <chrono>

#include "Bodies.cuh"
#include "Simulation.cuh"
#include "BoxBuilder.cuh"
#include "Engine.cuh"


#include <cuda.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <device_functions.h>
#include <cuda_runtime_api.h>





//const int THREADS_PER_MONITORBLOCK = 256;
//const int N_MONITORBLOCKS_PER_STEP = 256;

class Analyzer {
public:
	Analyzer() {}


	void analyzeEnergy(Simulation* simulation); // Prints a file of doubles: [step, molecule, atom, coordinate_dim]

	Float3* analyzeSolvateEnergy(Simulation* simulation, int n_steps);
	Float3* analyzeCompoundEnergy(Simulation* simulation, int n_steps);


private:
	void printEnergies(Float3* energy_data, int n_steps);
	Engine engine;
};
