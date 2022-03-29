#include "ForcefieldTypes.h"
#include "Filehandling.h"
#include "ForcefieldMerger.h"


#ifdef __linux__

	string sim_path = "/home/lima/Desktop/LIMA/Simulation";
#else
	string sim_path = "C:\\PROJECTS\\QUANTOM\\Simulation";
#endif

string mol_path = FileHelpers::pathJoin(sim_path, "Molecule");
string forcefield_path = FileHelpers::pathJoin(sim_path, "Forcefield");




vector<NB_Atomtype> makeFilteredNonbondedFF(Map* map) {
//	vector<vector<string>> simconf_rows = readFile("C:\\PROJECTS\\Quantom\\Molecules\\t4lys_full\\conf.gro");
	vector<vector<string>> simconf_rows = Reader::readFile(FileHelpers::pathJoin(mol_path, "conf.gro"));
	vector<string> simconf = FTHelpers::parseConf(simconf_rows);


//	vector<vector<string>> ffnonbonded_rows = readFile("C:\\PROJECTS\\Quantom\\charmm36-mar2019.ff\\ffnonbonded.itp");
	//vector<vector<string>> ffnonbonded_rows = readFile(sim_path+"/Forcefield/ffnonbonded.itp");
	vector<vector<string>> ffnonbonded_rows = Reader::readFile(FileHelpers::pathJoin(forcefield_path, "ffnonbonded.itp"));
	vector<NB_Atomtype> ffnonbonded = NB_Atomtype::parseNonbonded(ffnonbonded_rows);

	return NB_Atomtype::filterUnusedTypes(ffnonbonded, simconf, map);
}

vector<Atom> makeTopologyAtoms(vector<vector<string>> topology_rows, vector<NB_Atomtype>* ff_nonbonded_active, Map* map) {
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


void prepareForcefieldForSimulation() {
	// 95% of runtime is spent matching topologies to forcefield. 
	// To speedup, split the topology_x_types and run multithreaded?
	
	Map map;
	vector<NB_Atomtype> ff_nonbonded_active = makeFilteredNonbondedFF(&map);		// The map is made here, so this function must come first



	// These two vectors contain information about all types, but are parsed individually. This is very slightly slower, but MUCH more readable!
//	vector<vector<string>> ffbonded_rows = readFile("C:\\PROJECTS\\Quantom\\charmm36-mar2019.ff\\ffbonded.itp");
	//vector<vector<string>> ffbonded_rows = readFile(sim_path + "/Forcefield/ffbonded.itp");
	vector<vector<string>> ffbonded_rows = Reader::readFile(FileHelpers::pathJoin(forcefield_path, "ffbonded.itp"));
	//	vector<vector<string>> topology_rows = readFile("C:\\PROJECTS\\Quantom\\Molecules\\t4lys_full\\topol.top");
		//vector<vector<string>> topology_rows = readFile(sim_path + "/Molecule/topol.top");
	vector<vector<string>> topology_rows = Reader::readFile(FileHelpers::pathJoin(mol_path, "topol.top"));


	vector<Atom> atoms = makeTopologyAtoms(topology_rows, &ff_nonbonded_active, &map);



	vector<Bondtype> topology_bonds = makeTopologyBonds(&ffbonded_rows, &topology_rows, &atoms);

	vector<Angletype> topology_angles = makeTopologyAngles(&ffbonded_rows, &topology_rows, &atoms);

	vector<Dihedraltype> topology_dihedrals = makeTopologyDihedrals(ffbonded_rows, topology_rows, atoms);



	Printer::printForcefieldSummary(
		FileHelpers::pathJoin(forcefield_path, "ForcefieldSummary.txt"),
		ff_nonbonded_active, &map
	);


	Printer::printForcefield(

		FileHelpers::pathJoin(forcefield_path, "Forcefield.txt"),
		atoms,
		topology_bonds,
		topology_angles,
		topology_dihedrals
	);
}




int main(int argc, char* argv[]) {
	//prepareForcefieldForSimulation();


	vector<string> files;
	files.push_back("C:\\PROJECTS\\Quantom\\charmm36-mar2019.ff\\toppar_c36_jul18\\par_all35_ethers.prm");
	files.push_back("C:\\PROJECTS\\Quantom\\charmm36-mar2019.ff\\toppar_c36_jul18\\par_all36_carb.prm");
	files.push_back("C:\\PROJECTS\\Quantom\\charmm36-mar2019.ff\\toppar_c36_jul18\\par_all36_lipid.prm");
	files.push_back("C:\\PROJECTS\\Quantom\\charmm36-mar2019.ff\\toppar_c36_jul18\\par_all36_na.prm");
	files.push_back("C:\\PROJECTS\\Quantom\\charmm36-mar2019.ff\\toppar_c36_jul18\\par_all36m_prot.prm");
	
	

	files.push_back("C:\\PROJECTS\\Quantom\\charmm36-mar2019.ff\\toppar_c36_jul18\\par_all36_cgenff.prm");	// CHARMM wants this to be read last for some reason??


	ForcefieldMerger FM;
	FM.mergeForcefields(files);
	Printer::printFFNonbonded(FileHelpers::pathJoin(forcefield_path, "LIMA_ffnonbonded.txt"), FM.ff_nonbonded);
	Printer::printFFBonded(FileHelpers::pathJoin(forcefield_path, "LIMA_ffbonded.txt"), FM.ff_bondtypes, FM.ff_angletypes, FM.ff_dihedraltypes);
	
	return 0;
}





