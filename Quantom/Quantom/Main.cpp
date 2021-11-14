#include <iostream>
#include "Environment.h"
#include <math.h>


#include "Engine.cuh"

struct Test {
	int arr[4];
};
/*
int main() {
	Engine engine;

	return engine.testFunction();
}
*/


int main() {

	std::printf("Program starting...\n");
	Environment Env;

	Env.run();

	return 0;
}
