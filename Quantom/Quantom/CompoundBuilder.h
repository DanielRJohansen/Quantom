#pragma once


#include "Bodies.cuh"
#include <string.h>
#include <fstream>



using namespace std;


enum LineType { atom, pairbond, anglebond, torsionbond };

struct ParsedLine {
	LineType linetype;
	int pdb_index;
	string id;
	string monomer;
	Float3 position;
	string atom;
};





using namespace std;
class CompoundBuilder
{
public:
	CompoundBuilder() {}
	Compound buildMolecule(string pdb_path, string itp_path);




private:
	void loadParticles(Compound* compound, vector<vector<string>>* pdb_data, int max_monomer_id = INT32_MAX);
	void loadTopology(Compound* compound, vector<vector<string>>* itp_data, int* particle_id_map);


	enum TopologyMode { INACTIVE, BOND, ANGLE, DIHEDRAL };
	TopologyMode setMode(string entry);
	void addGeneric(Compound* compound, vector<string>* record, TopologyMode mode);
	void addBond(Compound* compound, vector<string>* record);
	void addAngle(Compound* compound, vector<string>* record);
	void addDihedral(Compound* compound, vector<string>* record);


	int* particle_id_map;


	ParsedLine parseLine(int line_index);
	ParsedLine parseAtom(string line);
	ParsedLine parseConnection(string line);
	vector<vector<string>> readFile(string path);
};

