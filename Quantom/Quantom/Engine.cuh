#pragma once

#include <iostream>

#include "Bodies.cuh"
#include "Simulation.cuh"



class Engine {
public:
	Engine(){}
	Engine(Simulation* simulation);

	void step();

private:
	Simulation* simulation;
	int fillBox();		// Returns # of bodies placed

	int initBlocks();	// returns # of blocks created

	void placeBody(SimBody* body);





	// Helper functions
	Int3 posToBlockIndex(Double3* pos) {
		return Int3(
			floor((pos->x + simulation->box_size / 2.f) / simulation->box_size * simulation->blocks_per_dim),
			floor((pos->y + simulation->box_size / 2.f) / simulation->box_size * simulation->blocks_per_dim),
			floor((pos->z + simulation->box_size / 2.f) / simulation->box_size * simulation->blocks_per_dim)
			);
	}


	bool finished = false;
};