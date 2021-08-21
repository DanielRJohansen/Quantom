#pragma once

//#include "QuantomTypes.cuh"
#include "Bodies.cuh"
#include "Display.h"
#include "Interface.h"


class Environment
{
public:
	Environment();





private:
	Display* display;
	Interface* interface;

	void run();

};

