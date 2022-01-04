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
	/*
	uint8_t* a;
	a = new uint8_t[4];

	void* b = a;
	//*(uint32_t*)b = 0xFF000100;
	//*(uint32_t*)b ^= *(uint32_t*)b;
	//*(uint32_t*)a = 0xFF000000;
	uint32_t c = 0xF0000000;
	
	*a = *(uint8_t*)c;
	for (int i = 0; i < 4; i++)
		printf("%d\n", a[i]);

	
	exit(0);
	*/


	std::printf("Program starting...\n");
	Environment Env;

	Env.run();
	//Env.makeVirtualTrajectory("D:\\Quantom\\trajectory.csv", "D:\\Quantom\\waterforce.csv");
	//Env.renderTrajectory("D:\\Quantom\\virtrj.csv");
	//Env.renderTrajectory("D:\\Quantom\\trajectory.csv");
	return 0;
}
