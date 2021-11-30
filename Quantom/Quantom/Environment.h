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






class Environment
{
public:
	Environment();

	void run();
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



	// These should be in interface maybe?
	void printOut(double* data, int n_steps);
	//void printDataBuffer(Box* box);
	void printTrajectory(Simulation* simulation);
	void printWaterforce(Simulation* simulation);
};

