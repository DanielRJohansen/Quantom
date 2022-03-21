#pragma once

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>

///////////////////////////////// READ HERE FIRST /////////////////////////////////
// ffbonded.itp and ffnonbonden.itp has different atomtypes. 
// The purpose if this bit of code is to compress all forcefield parameters so they can fit on a GPU more easily
//

using namespace std;







struct Map {
	struct Mapping {
		Mapping () {}
		Mapping (string l, string r) : left(l), right(r) {}
		string left, right;
	};



	Map() {
		mappings = new Mapping[10000];
	}
	~Map() {
		delete[] mappings;
	}

	Mapping* mappings;
	int n_mappings = 0; 

	void addMap(string l, string r) {
		mappings[n_mappings++] = Mapping(l, r);
	}
	string mapRight(string query) {
		for (int i = 0; i < n_mappings; i++) {
			if (mappings[i].left == query) {
				return mappings[i].right;
			}

		}
		return query;
		//cout << "No map for query" << query << endl;
		//exit(0);
	}
};


struct FF_nonbonded {
	FF_nonbonded() {}
	FF_nonbonded(string t) : type(t){}		// This is for loading form the conf file, for easy comparisons
	FF_nonbonded(string t, int a, double m, double s, double e) : type(t), atnum(a), mass(m), sigma(s), epsilon(e) {}
	// Official parameters
	string type;
	int atnum = -1;		// atnum given by CHARMM
	int atnum_local;	// atnum specific to simulation
	double mass, sigma, epsilon;

	// LIMA parameters
	bool is_present_in_simulation = false;
	int simulation_specific_id;



	//bool operator==(const FF_nonbonded a) { return (type == a.type); }
	//static bool sameType(FF_nonbonded a, FF_nonbonded b) { return (a.type == b.type); }


	// Parser functions
	enum STATE { INACTIVE, ATOMTYPES, PAIRTYPES };
	static STATE setState(string s, STATE current_state) {
		if (s == "atomtypes")
			return ATOMTYPES;
		if (s == "pairtypes")
			return PAIRTYPES;
		return current_state;
	}
	static vector<FF_nonbonded> parseNonbonded(vector<vector<string>> rows) {		// KNOWN ERROR HERE. SOMETIMES SPACES IN FILES ARE TABS, AND IT DOESN'T RECOGNISE A ROW AS A PROPER ENTRY!!
		STATE current_state = INACTIVE;

		vector<FF_nonbonded> records;

		for (vector<string> row : rows) {

			if (row.size() == 3) {
				current_state = setState(row[1], current_state);
				continue;
			}


		START:

			switch (current_state)
			{
			case INACTIVE:
				break;
			case ATOMTYPES:
				if (row.size() != 7) {
					continue;													// TODO: dont just run from your problems, man...

				}
				records.push_back(FF_nonbonded(row[0], stoi(row[1]), stod(row[2]), stod(row[5]), stod(row[6])));
				break;
			case PAIRTYPES:
				break;
			default:
				break;
			}
		}
		printf("%d atom types read from file\n", records.size());
		return records;
	}
	static FF_nonbonded findRecord(vector<FF_nonbonded>* records, string type) {
		for (FF_nonbonded record : *records) {
			if (record.type == type) {
				return record;
			}
		}
		return FF_nonbonded();
	}
	static bool typeAlreadyPresent(vector<FF_nonbonded>* records, string type) {
		return (findRecord(records, type).atnum != -1);
	}
	static vector<FF_nonbonded> filterUnusedTypes(vector<FF_nonbonded> records, vector<string> active_types, Map* map) {
		vector<FF_nonbonded> filtered_list;
		//string* mappings = new string[10000];
		for (string type : active_types) {
			bool found = false;
			string alias = type;
			TRY_AGAIN:



			FF_nonbonded record = findRecord(&records, alias);

			if (record.atnum != -1) {
				if (!typeAlreadyPresent(&filtered_list, alias)) {
					//cout << "Adding alias " << alias << endl;
					record.atnum_local = filtered_list.size();
					filtered_list.push_back(record);
				}
				if (type != alias) {
					map->addMap(type, alias);
					//mappings[*m_cnt * 2] = type;
					//mappings[*m_cnt * 2 + 1] = alias;
					//(*m_cnt)++;
				}
						
									
			}
			else {
				//cout << "Atom type " << alias << " not found" << endl;

				if (alias.length() > 0) {
					alias.pop_back();
					goto TRY_AGAIN;
				}
				else {
					printf("Failed");
					exit(404);
				}
			}
		}
		printf("Filtered ff down to %d entries (%d bytes)\n", filtered_list.size(), sizeof(float) * 3 * filtered_list.size());
		printf("Aliases (%d): \n", map->n_mappings);
		for (int i = 0; i < map->n_mappings; i++) {
			//cout << mappings[i*2] << "    " << mappings[i*2+1] << endl;
		}
		return filtered_list;
	}
	
};

struct FF_bondtype {
	FF_bondtype() {}
	FF_bondtype(string t1, string t2) : type1(t1), type2(t2) {
		sort();
	}
	FF_bondtype(string t1, string t2, float b0, float kb) : type1(t1), type2(t2), b0(b0), kb(kb) {
		sort();
	}
	FF_bondtype(int id1, int id2) : id1(id1), id2(id2) {
		//convertToZeroIndexed();
	}
	
	string type1, type2;
	int id1, id2;			// bonds from .top only has these values! 
	float b0;
	float kb;

	void convertToZeroIndexed() {
		id1--;
		id2--;
	}
	void sort() {
		int ptr = 0; 
		while (type1.length() > ptr && type2.length() > ptr) {
			if ((int)type1[ptr] > (int)type2[ptr]) {
				swap(type1, type2);
				return;
			}
			ptr++;
		}
		if (type1.length() > type2.length()) {
			swap(type1, type2);
		}
	}

	enum STATE { INACTIVE, BONDTYPES, ANGLETYPES, DIHEDRALTYPES };
	static STATE setState(string s, STATE current_state) {
		if (s == "bondtypes" || s == "bonds")
			return BONDTYPES;
		if (s == "angletypes")
			return ANGLETYPES;
		if (s == "dihedraltypes")
			return DIHEDRALTYPES;
		return current_state;
	}
	static vector<FF_bondtype> parseFFBondtypes(vector<vector<string>> rows) {
		STATE current_state = INACTIVE;

		vector<FF_bondtype> records;

		for (vector<string> row : rows) {

			if (row.size() == 3) {
				current_state = setState(row[1], current_state);
				continue;
			}



			switch (current_state)
			{
			case INACTIVE:
				break;
			case BONDTYPES:
				if (row.size() != 5) {
					continue;													// TODO: dont just run from your problems, man...
				}
				records.push_back(FF_bondtype(row[0], row[1], stof(row[3]), stof(row[4])));
				break;
			case ANGLETYPES:
				break;
			case DIHEDRALTYPES:
				break;
			default:
				break;
			}
		}
		printf("%d single bonds read from file\n", records.size());
		return records;
	}

	static vector<FF_bondtype> parseTopolBondtypes(vector<vector<string>> rows) {
		STATE current_state = INACTIVE;
		vector<FF_bondtype> records;

		for (vector<string> row : rows) {
			if (row.size() == 3) {
				if (row[0][0] == '[') {
					current_state = setState(row[1], current_state);
					continue;
				}				
			}



			switch (current_state)
			{
			case BONDTYPES:
				records.push_back(FF_bondtype(stoi(row[0]), stoi(row[1])));
				break;
			default:
				break;
			}

			if (current_state == ANGLETYPES)
				break;
		}
		printf("%d bonds found in topology file\n", records.size());
		return records;	
	}

	static void assignTypesFromAtomIDs(vector<FF_bondtype>* topol_bonds, vector<FF_bondtype> FF_bondtypes, vector<Atom> atoms) {
		for (int i = 0; i < topol_bonds->size(); i++) {
			FF_bondtype* bond
		}
	}
};


// This is for bonded atoms!!!!!!!!!!!
struct Atom {
	Atom() {}
	Atom(int id, string type_b, string type_nb) : id(id), atomtype_bond(type_b), atomtype(type_nb) {
		//convertToZeroindexed();
	}
	int id;										// Come from topol.top file
	string atomtype;
	string atomtype_bond;	
	int atomtype_id = -1;				// Asigned later
	//float charge;


	void convertToZeroindexed() { id--; }

	enum STATE { INACTIVE, ATOMS, FINISHED};
	static STATE setState(string s, STATE current_state) {
		if (s == "atoms")
			return ATOMS;
		if (s == "bonds")
			return FINISHED;
		return current_state;
	}

	static vector<Atom> parseTopolAtoms(vector<vector<string>> rows) {
		STATE current_state = INACTIVE;

		vector<Atom> records;

		for (vector<string> row : rows) {

			if (row.size() == 3) {
				current_state = setState(row[1], current_state);
				continue;
			}



			switch (current_state)
			{
			case INACTIVE:
				break;
			case ATOMS:
				records.push_back(Atom(stoi(row[0]), row[1], row[4]));
				break;
			default:
				break;
			}

			if (current_state == FINISHED)
				break;
		}
		printf("%d atoms found in topology file\n", records.size());
		return records;
	}

	//void assignAtomtypeID()
	static void assignAtomtypeIDs(vector<Atom>* atoms, vector<FF_nonbonded>* forcefield, Map* map) {
		

		for (int i = 0; i < atoms->size(); i++) {
			Atom* atom = &((*atoms).at(i));
			string alias = map->mapRight(atom->atomtype);

			for (FF_nonbonded force_parameters : *forcefield) {
				if (force_parameters.type == alias) {
					atom->atomtype_id = force_parameters.atnum_local;
				}
			}
		}

		/*
		for (Atom atom : *atoms) {
			string alias = map->mapRight(atom.atomtype);

			for (FF_nonbonded force_parameters : *forcefield) {
				if (force_parameters.type == alias) {
					atom.atomtype_id = force_parameters.simulation_specific_id;
				}
			}
		}
		*/

		for (Atom atom : *atoms) {
			if (atom.atomtype_id == -1) {
				cout << atom.atomtype << "   ";
				printf("atom id %d\n", atom.atomtype_id);
				exit(0);
			}
		}
	}
};





static vector<string> parseConf(vector<vector<string>> rows) {		// KNOWN ERROR HERE. SOMETIMES SPACES IN FILES ARE TABS, AND IT DOESN'T RECOGNISE A ROW AS A PROPER ENTRY!!
	vector<string> atom_types;

	for (vector<string> row : rows) {
		if (row.size() != 6) {
			continue;
		}

		bool type_already_found = false;
		for (string atom_type : atom_types) {
			if (atom_type == row[1]) {
				type_already_found = true;
				break;
			}
		}
		if (!type_already_found)
			atom_types.push_back(row[1]);
	}

	printf("%d atom types in conf\n", atom_types.size());

	return atom_types;
}



