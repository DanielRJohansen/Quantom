#pragma once

#include <iostream>

#include "Bodies.cuh"
#include "Simulation.h"



class Engine {
public:
	Engine(){}
	Engine(Simulation* simulation);

	void step();

private:
	Simulation* simulation;
	void fillBox();




	bool finished = false;
};