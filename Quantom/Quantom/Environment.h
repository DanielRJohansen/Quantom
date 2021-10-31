#pragma once

//#include "QuantomTypes.cuh"
#include "Bodies.cuh"
#include "Display.h"
#include "Interface.h"
#include "Engine.cuh"

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

	void printOut(float* data);

};

