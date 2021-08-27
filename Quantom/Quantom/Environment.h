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





private:
	Display* display;
	Interface* interface;
	Simulation simulation;
	Engine* engine;
	void run();

};

