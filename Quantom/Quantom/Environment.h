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




#include <stdlib.h>
#include <stdio.h>




#ifndef __linux__
#include <direct.h>
#endif










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
	void verifySimulationParameters();			// Constants before doing anything
	void verifyBox();							// Checks wheter the box will break
	//Display* display;
	DisplayV2* display;
	Interface* interface;
	Engine* engine;
	ForceFieldMaker* forcefieldmaker;
	Analyzer analyzer;
	CompoundBuilder* compoundbuilder;
	BoxBuilder boxbuilder;




	// These should be in interface maybe?
	template <typename T>
	void dumpToFile(T* data, int n_datapoints, string file_path);



#ifdef __linux__
	std::chrono::system_clock::time_point time0;
#else
	std::chrono::steady_clock::time_point time0;
#endif

};

