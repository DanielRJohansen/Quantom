#pragma once

#include "QuantomTypes.cuh"
#include <SFML/Graphics.hpp>
#include <iostream>

enum Atom{Oxygen};

__global__ struct Sphere {
	Sphere() {}
	Sphere(Double3 pos, double r, sf::Color c, double mass) : pos(pos), radius(r), color(c), mass(mass) {}
	Double3 pos;
	double radius;
	double mass;
	sf::Color color;
};

__global__ struct RenderBody {

};

__global__ struct Molecule {
	Molecule();
	int n_atoms;
	Sphere* atoms;
	Double3 CoM;
};


__global__ struct SimBody {
	//RenderBody* renderbody;



	Double3 pos;	//CoM
	Double3 vel;
	double mass;



	Double3 charge_unit_vector;
	double charge_magnitude;

};



__global__ struct Compound {	//Or molecule
	Compound() {};
	Compound(int n_bodies) : n_bodies(n_bodies) {
		bodies = new SimBody[n_bodies]; 
	}

	bool non_bonded = true;

	int n_bodies = 0;
	SimBody* bodies;
};