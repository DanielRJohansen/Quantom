#pragma once

//#include "QuantomTypes.cuh"
#include "Bodies.cuh"
#include "Display.h"
#include "Interface.h"
#include "Engine.cuh"

class Environment
{
public:
	Environment();

	void run();



	Simulation* simulation;


private:
	Display* display;
	Interface* interface;
	Engine* engine;

};

