#include "ForcefieldTypes.h"
#include "Filehandling.h"














int main(int argc, char* argv[]) {
	vector<vector<string>> ffnonbonded_rows = readFile("C:\\PROJECTS\\Quantom\\charmm36-mar2019.ff\\ffnonbonded.itp", 2000);
	vector<FF_nonbonded> ffnonbonded = FF_nonbonded::parseNonbonded(ffnonbonded_rows);


	vector<vector<string>> simconf_rows = readFile("C:\\PROJECTS\\Quantom\\Molecules\\t4lys_full\\conf.gro", 2000);
	vector<FF_nonbonded> simconf = FF_nonbonded::parseConf(simconf_rows);



	printf("Hello world!");

	return 0;
}





