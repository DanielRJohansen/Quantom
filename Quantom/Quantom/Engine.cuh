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
	void fillBox();
	void initBlocks();

	void placeBody(SimBody* body);





	// Helper functions
	Int3 posToBlockIndex(Double3* pos) {
		return Int3(
			(pos->x + simulation->box_size / 2) / simulation->box_size * simulation->blocks_per_dim,
			(pos->y + simulation->box_size / 2) / simulation->box_size * simulation->blocks_per_dim,
			(pos->z + simulation->box_size / 2) / simulation->box_size * simulation->blocks_per_dim
			);
	}


	bool finished = false;
};