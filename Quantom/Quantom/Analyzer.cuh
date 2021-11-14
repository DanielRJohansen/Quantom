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




class Analyzer {
public:
	Analyzer() {}


	void analyzeEnergy(Simulation* simulation); // Prints a file of floats: [step, molecule, atom, coordinate_dim]



private:
	void printEnergies(float* energy_data, int n_steps);
	Engine engine;
};
