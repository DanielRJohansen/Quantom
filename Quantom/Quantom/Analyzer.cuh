#pragma once

#include <iostream>
#include <chrono>
#include <cstring>
#include <string>
//#include "Bodies.cuh"
#include "Simulation.cuh"
//#include "BoxBuilder.cuh"
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

	struct AnalyzedPackage {
		AnalyzedPackage(Float3* e_ptr, int e_cnt, float* t_ptr, int t_cnt) : n_energy_values(e_cnt), n_temperature_values(t_cnt) {
			energy_data = e_ptr;
			temperature_data = t_ptr;
		}
		Float3* energy_data; // potE, kinE, totalE
		int n_energy_values;

		float* temperature_data;
		int n_temperature_values;
		~AnalyzedPackage() {
			delete[] energy_data;
			delete[] temperature_data;
		}
	};

	AnalyzedPackage analyzeEnergy(Simulation* simulation); // Prints a file of doubles: [step, molecule, atom, coordinate_dim]

	Float3* analyzeSolvateEnergy(Simulation* simulation, int n_steps);
	Float3* analyzeCompoundEnergy(Simulation* simulation, int n_steps);



private:
	void printEnergies(Float3* energy_data, int n_steps, Simulation* simulation);
	Engine engine;



	Float3* traj_buffer_device;
	double* potE_buffer_device;
};
