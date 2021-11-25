#pragma once

//#include "QuantomTypes.cuh"
#include "Bodies.cuh"
#include "Display.h"
#include "Interface.h"
#include "Engine.cuh"
#include "Analyzer.cuh"
#include "CompoundBuilder.h"


// For logging
#include <fstream>
#include <string>



class Environment
{
public:
	Environment();

	void run();



	Simulation* simulation;


private:
	bool verifySimulationParameters();
	Display* display;
	Interface* interface;
	Engine* engine;
	Analyzer analyzer;
	CompoundBuilder compoundbuilder;


	// These should be in interface maybe?
	void printOut(float* data);
	void printDataBuffer(Box* box);
	void printTrajectory(Simulation* simulation);
};

