#include <iostream>
#include "Environment.h"
#include <math.h>



int main() {

/*
	SimBody bodies[19];
	bodies[10] = SimBody();
	bodies[10].pos.print();
	//bodies[10].molecule_type = 10;
	printf("%d\n",  bodies[10].molecule_type);
	printf("%f", bodies[10].charge_magnitude);
	exit(1);
	*/

	std::printf("Program starting...\n");
	Environment Env;

	Env.run();

	return 1;
}
