#include "ForcefieldTypes.h"
#include "Filehandling.h"



#ifdef __linux__
	string sim_path = "/home/lima/Desktop/LIMA/Simulation";
#else
	string sim_path = "C:\\PROJECTS\\QUANTOM\\Molecules";
#endif






vector<FF_nonbonded> makeFilteredNonbondedFF(Map* map) {
//	vector<vector<string>> simconf_rows = readFile("C:\\PROJECTS\\Quantom\\Molecules\\t4lys_full\\conf.gro");
	vector<vector<string>> simconf_rows = readFile(sim_path + "/Molecule/conf.gro");
	vector<string> simconf = parseConf(simconf_rows);


//	vector<vector<string>> ffnonbonded_rows = readFile("C:\\PROJECTS\\Quantom\\charmm36-mar2019.ff\\ffnonbonded.itp");
	vector<vector<string>> ffnonbonded_rows = readFile(sim_path+"/Forcefield/ffnonbonded.itp");
	vector<FF_nonbonded> ffnonbonded = FF_nonbonded::parseNonbonded(ffnonbonded_rows);

	return FF_nonbonded::filterUnusedTypes(ffnonbonded, simconf, map);
}

vector<Atom> makeTopologyAtoms(vector<vector<string>> topology_rows, vector<FF_nonbonded>* ff_nonbonded_active, Map* map) {
	vector<Atom> atoms = Atom::parseTopolAtoms(topology_rows);
	Atom::assignAtomtypeIDs(&atoms, ff_nonbonded_active, map);
	return atoms;
}

vector<Bondtype> makeTopologyBonds(vector<vector<string>>* ffbonded_rows, vector<vector<string>>* topology_rows, vector<Atom>* atoms) {
	vector<Bondtype> ffbondtypes = Bondtype::parseFFBondtypes(*ffbonded_rows);
	vector<Bondtype> topology_bonds = Bondtype::parseTopolBondtypes(*topology_rows);

	Bondtype::assignTypesFromAtomIDs(&topology_bonds, *atoms);
	Bondtype::assignFFParametersFromBondtypes(&topology_bonds, &ffbondtypes);
	
	return topology_bonds;
}

vector<Angletype> makeTopologyAngles(vector<vector<string>>* ffbonded_rows, vector<vector<string>>* topology_rows, vector<Atom>* atoms) {
	vector<Angletype> ffangletypes = Angletype::parseFFAngletypes(*ffbonded_rows);
	vector<Angletype> topology_angles = Angletype::parseTopolAngletypes(*topology_rows);

	Angletype::assignTypesFromAtomIDs(&topology_angles, *atoms);
	Angletype::assignFFParametersFromAngletypes(&topology_angles, &ffangletypes);

	return topology_angles;
}

vector<Dihedraltype> makeTopologyDihedrals(vector<vector<string>> ffbonded_rows, vector<vector<string>> topology_rows, vector<Atom> atoms) {
	vector<Dihedraltype> forcefield = Dihedraltype::parseFFDihedraltypes(ffbonded_rows);
	vector<Dihedraltype> topology_dihedrals = Dihedraltype::parseTopolDihedraltypes(topology_rows);

	Dihedraltype::assignTypesFromAtomIDs(&topology_dihedrals, atoms);
	Dihedraltype::assignFFParametersFromDihedraltypes(&topology_dihedrals, &forcefield);

	return topology_dihedrals;
}

int main(int argc, char* argv[]) {
	// 95% of runtime is spent matching topologies to forcefield. 
	// To speedup, split the topology_x_types and run multithreaded?
	//


	Map map;
	vector<FF_nonbonded> ff_nonbonded_active = makeFilteredNonbondedFF(&map);		// The map is made here, so this function must come first
	


	// These two vectors contain information about all types, but are parsed individually. This is very slightly slower, but MUCH more readable!
//	vector<vector<string>> ffbonded_rows = readFile("C:\\PROJECTS\\Quantom\\charmm36-mar2019.ff\\ffbonded.itp");
	vector<vector<string>> ffbonded_rows = readFile(sim_path + "/Forcefield/ffbonded.itp");
//	vector<vector<string>> topology_rows = readFile("C:\\PROJECTS\\Quantom\\Molecules\\t4lys_full\\topol.top");
	vector<vector<string>> topology_rows = readFile(sim_path + "/Molecule/topol.top");


	vector<Atom> atoms = makeTopologyAtoms(topology_rows, &ff_nonbonded_active, &map);


	
	vector<Bondtype> topology_bonds = makeTopologyBonds(&ffbonded_rows, &topology_rows, &atoms);

	vector<Angletype> topology_angles = makeTopologyAngles(&ffbonded_rows, &topology_rows, &atoms);

	vector<Dihedraltype> topology_dihedrals = makeTopologyDihedrals(ffbonded_rows, topology_rows, atoms);
	














	printForcefieldSummary(
//		(string)"C:\\Users\\Daniel\\git_repo\\Quantom\\" + (string)"ForcefieldSummary.txt",
		sim_path + "/Forcefield/ForcefieldSummary.txt",
		ff_nonbonded_active, &map
	);
	
	
	printForcefield(
//		"C:\\Users\\Daniel\\git_repo\\Quantom\\" + (string)"Forcefield.txt",
		sim_path + "/Forcefield/Forcefield.txt",
		atoms,
		topology_bonds,
		topology_angles,
		topology_dihedrals
	);
	

	printf("Hello world!");

	return 0;
}





