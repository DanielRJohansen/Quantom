#include <iostream>
#include "Environment.h"
#include <math.h>


#include "Engine.cuh"

struct Test {

	int arr[4];
};



int main() {


	std::printf("Program starting...\n");
	Environment Env;

	Env.run();
	//Env.makeVirtualTrajectory("D:\\Quantom\\trajectory.csv", "D:\\Quantom\\waterforce.csv");
	//Env.renderTrajectory("D:\\Quantom\\virtrj.csv");
	//Env.renderTrajectory("D:\\Quantom\\trajectory.csv");
	return 0;
}
