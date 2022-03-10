#pragma once


#include "Bodies.cuh"
#include <string.h>
#include <fstream>
#include <vector>

#include "Forcefield.cuh"


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
	Compound buildMolecule(string pdb_path, string itp_path, int max_residue_id=INT16_MAX);




private:
	ForceFieldMaker FFM;


	struct Record_ATOM;
	void loadParticles(Compound* compound, vector<Record_ATOM>* pdb_data, int max_monomer_id = INT32_MAX, bool ignore_protons =false);
	void loadTopology(Compound* compound, vector<vector<string>>* itp_data, int* particle_id_map);
	void calcParticleSphere(Compound* compound);


	enum TopologyMode { INACTIVE, BOND, ANGLE, DIHEDRAL };
	TopologyMode setMode(string entry);
	void addGeneric(Compound* compound, vector<string>* record, TopologyMode mode);
	void addBond(Compound* compound, vector<string>* record);
	void addAngle(Compound* compound, vector<string>* record);
	void addDihedral(Compound* compound, vector<string>* record);


	Float3 calcCOM(Compound* compound);

	int* particle_id_map;


	ParsedLine parseLine(int line_index);
	ParsedLine parseAtom(string line);
	ParsedLine parseConnection(string line);

	vector<vector<string>> parseTOP(string path);
	vector<Record_ATOM> parsePDB(string path);
	vector<Record_ATOM> parseGRO(string path);


	struct ResidueComboId;

	ResidueComboId parseResidueID(string s);
	bool isAsciiNumber(char c);




	struct ResidueComboId {
		ResidueComboId(){}
		ResidueComboId(int id, string name) : id(id), name(name) { valid = true; }
		bool valid = false;
		int id;
		string name;
	};





	struct Record_ATOM {	//	https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
		Record_ATOM(int a, string b, char c, string d, char e, int f, char g, Float3 h) {
			atom_serial_number = a;
			atom_name = b;
			alternate_location_indicator = c;
			residue_name = d;
			chain_identifier = e;
			residue_seq_number = f;
			code_for_insertions_of_residues = g;
			coordinate = h;
		}

		int atom_serial_number;
		string atom_name;
		char alternate_location_indicator;
		string residue_name;
		char chain_identifier;
		int residue_seq_number;
		char code_for_insertions_of_residues;
		Float3 coordinate;						// [nm]
	};




};





