#include "ForcefieldTypes.h"
#include "Filehandling.h"














int main(int argc, char* argv[]) {
	vector<vector<string>> simconf_rows = readFile("C:\\PROJECTS\\Quantom\\Molecules\\t4lys_full\\conf.gro");
	vector<string> simconf = parseConf(simconf_rows);
	
	
	vector<vector<string>> ffnonbonded_rows = readFile("C:\\PROJECTS\\Quantom\\charmm36-mar2019.ff\\ffnonbonded.itp");
	vector<FF_nonbonded> ffnonbonded = FF_nonbonded::parseNonbonded(ffnonbonded_rows);


	Map map;
	vector<FF_nonbonded> ff_nonbonded_active = FF_nonbonded::filterUnusedTypes(ffnonbonded, simconf, &map);
	





	vector<vector<string>> ffbonded_rows = readFile("C:\\PROJECTS\\Quantom\\charmm36-mar2019.ff\\ffbonded.itp");
	vector<FF_bondtype> ffbondtypes = FF_bondtype::parseFFBondtypes(ffbonded_rows);

	

	// Now for atoms and bonds present in topology
	vector<vector<string>> topology_rows = readFile("C:\\PROJECTS\\Quantom\\Molecules\\t4lys_full\\topol.top");
	vector<Atom> atoms = Atom::parseTopolAtoms(topology_rows);
	vector<FF_bondtype> topology_bonds = FF_bondtype::parseTopolBondtypes(topology_rows);


	Atom::assignAtomtypeIDs(&atoms, &ff_nonbonded_active, &map);


	FF_bondtype::assignTypesFromAtomIDs(&topology_bonds, atoms);
	FF_bondtype::assignFFParametersFromBondtypes(&topology_bonds, &ffbondtypes);



	printForcefieldSummary(
		(string)"C:\\Users\\Daniel\\git_repo\\Quantom\\" + (string)"ForcefieldSummary.txt",
		ff_nonbonded_active, &map
	);
	printForcefield(
		"C:\\Users\\Daniel\\git_repo\\Quantom\\" + (string)"Forcefield.txt",
		atoms,
		topology_bonds
	);
	

	printf("Hello world!");

	return 0;
}





