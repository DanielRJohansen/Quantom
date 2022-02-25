#pragma once

//#include "QuantomTypes.cuh"
#include "Bodies.cuh"
//#include "Display.h"
#include "DisplayV2.h"
#include "Interface.h"
#include "Engine.cuh"
#include "Analyzer.cuh"
#include "CompoundBuilder.h"
#include "VirtualPathMaker.h"
#include "BoxBuilder.cuh"


// For logging
#include <fstream>
#include <string>
#include <assert.h>  





class Environment
{
public:
	Environment();

	void run();
	void postRunEvents();
	void handleStatus(Simulation* simulation);
	void handleDisplay(Simulation* simulation);
	bool handleTermination(Simulation* simulation);


	void renderTrajectory(string trj_path);
	void makeVirtualTrajectory(string trj_path, string waterforce_path);

	Simulation* simulation;


private:
	bool verifySimulationParameters();
	//Display* display;
	DisplayV2* display;
	Interface* interface;
	Engine* engine;
	Analyzer analyzer;
	CompoundBuilder compoundbuilder;
	BoxBuilder boxbuilder;




	// These should be in interface maybe?
	//void printOut(double* data, int n_steps);
	template <typename T>
	void dumpToFile(T* data, int n_datapoints, string file_path);
	//void dumpToFile(double* data, Int3 dim, string file_name);
	//void dumpToFile(Float3* data, Int3 dim, string file_name);
	//void printDataBuffer(Box* box);
	void printTrajectory(Simulation* simulation);
	void printWaterforce(Simulation* simulation);

	void printFloat3Matrix(Float3* data_matrix, Int3 dim, string filename);






	std::chrono::steady_clock::time_point time0;

};

