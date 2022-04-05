#pragma once


#include "Bodies.cuh"
#include <string.h>
#include <fstream>
#include <vector>

#include "Constants.cuh"
#include "Forcefield.cuh"


using namespace std;


enum LineType { atom, singlebond, anglebond, torsionbond };

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
	struct IDMap {				// Delete for particle refs!
		IDMap() {}	
		IDMap(int global_id, int compound_id, int local_id) : global_id(global_id), compound_id(compound_id), local_id(local_id) {}
		int global_id = -1, compound_id = -1, local_id = -1;	// local_id is local to the compound
	};
public:
	CompoundBuilder() {}
	CompoundBuilder(ForceFieldMaker* ffm) { FFM = ffm; }
	Molecule buildMolecule(string pdb_path, string itp_path, int max_residue_id=INT16_MAX, int min_residue_id=0);




private:
	ForceFieldMaker* FFM;
	//IDMap* particle_id_maps;
	ParticleRef* particle_id_maps;
	CompoundBridgeBundle* compound_bridge_bundle;


	struct Record_ATOM;
	void loadParticles(Molecule* molecule, vector<Record_ATOM>* pdb_data, int max_monomer_id = INT32_MAX, int min_residue_id=0, bool ignore_protons =false);
	void loadTopology(Molecule* molecule, vector<vector<string>>* itp_data);
	void calcParticleSphere(Compound* compound);


	enum TopologyMode { INACTIVE, BOND, ANGLE, DIHEDRAL };
	TopologyMode setMode(string entry);
	void loadMaps(ParticleRef* maps, vector<string>* record, int n);
	void addGeneric(Molecule* molecule, vector<string>* record, TopologyMode mode);
	void addBond(Molecule* molecule, ParticleRef* maps, vector<string>* record);
	void addAngle(Molecule* molecule, ParticleRef* maps, vector<string>* record);
	void addDihedral(Molecule* molecule, ParticleRef* maps, vector<string>* record);
	void distributeLJIgnores(Molecule* molecule, ParticleRef* maps, int n);



	ParsedLine parseLine(int line_index);
	ParsedLine parseAtom(string line);
	ParsedLine parseConnection(string line);

	vector<vector<string>> parseTOP(string path);
	vector<Record_ATOM> parsePDB(string path);
	vector<Record_ATOM> parseGRO(string path);


	struct ResidueComboId;

	ResidueComboId parseResidueID(string s);
	inline bool isAsciiNumber(char c) { return (c > 47 && c < 58); }
	void countElements(Molecule* molecule);





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





