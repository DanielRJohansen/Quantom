#pragma once

//#include "QuantomTypes.cuh"
#include "Bodies.cuh"
#include "Display.h"
#include "Interface.h"
#include "Engine.cuh"
#include "Analyzer.cuh"
#include "CompoundBuilder.h"
#include "VirtualPathMaker.h"

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
	Display* display;
	Interface* interface;
	Engine* engine;
	Analyzer analyzer;
	CompoundBuilder compoundbuilder;
	BoxBuilder boxbuilder;




	// These should be in interface maybe?
	void printOut(double* data, int n_steps);
	//void printDataBuffer(Box* box);
	void printTrajectory(Simulation* simulation);
	void printWaterforce(Simulation* simulation);








	std::chrono::steady_clock::time_point time0;

};

