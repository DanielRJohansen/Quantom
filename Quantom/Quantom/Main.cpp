#include <iostream>
#include "Environment.h"
#include <math.h>


#include "Engine.cuh"

struct Test {

	int arr[4];
	//const static int size = 14;
};



int main() {


	//Environment::prepFF("conf_test.gro", "topol_test.top");

	string conf_filename = "conf.gro";
	string topol_filename = "topol.top";
	//Environment::prepFF("conf.gro", "topol.top");
	Environment::prepFF(conf_filename, topol_filename);
	std::printf("Program starting...\n");

	Environment Env(conf_filename, topol_filename);
	Env.run();

	//Env.makeVirtualTrajectory("D:\\Quantom\\trajectory.csv", "D:\\Quantom\\waterforce.csv");
	//Env.renderTrajectory("D:\\Quantom\\virtrj.csv");
	//Env.renderTrajectory("D:\\Quantom\\trajectory.csv");
	return 0;
}
