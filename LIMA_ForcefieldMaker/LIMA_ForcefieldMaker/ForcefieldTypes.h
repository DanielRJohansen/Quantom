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



bool charIsNumber(char c) {
	return ((int)c > 47 && (int)c < 58);
}
bool charIsNumberAbove1(char c) {
	return ((int)c > 49 && (int)c < 58);
}
float calcLikeness(string a, string b) {
	float likeness = 0;
	float point_scale = 1.f / max(a.length(), b.length());

	for (int i = 0; i < min(a.length(), b.length()); i++) {
		if (a[i] == b[i])
			likeness += point_scale;
		else
			break;
	}
	return likeness;
}

bool isSorted(string* leftmost, string* rightmost) {
	int ptr = 0;
	while (leftmost->length() > ptr && rightmost->length() > ptr) {
		if ((int)(*leftmost)[ptr] < (int)(*rightmost)[ptr]) {
			return true;
		}
		else if ((int)(*leftmost)[ptr] > (int)(*rightmost)[ptr]) {
			return false;
		}
		ptr++;
	}
	if (leftmost->length() > rightmost->length()) {
		return false;
	}
}

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

	enum STATE { INACTIVE, ATOMS, FINISHED };
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



struct Bondtype {
	Bondtype() {}
	Bondtype(string t1, string t2) : type1(t1), type2(t2) {
		sort();
	}
	Bondtype(string t1, string t2, float b0, float kb) : type1(t1), type2(t2), b0(b0), kb(kb) {
		sort();
	}
	Bondtype(int id1, int id2) : id1(id1), id2(id2) {
		//convertToZeroIndexed();
	}
	
	string type1, type2;
	int id1, id2;			// bonds from .top only has these values! 
	float b0;
	float kb;

	
	void sort() {
		int ptr = 0; 
		//cout << "Sorting " << type1 << "    " << type2 << endl;
		while (type1.length() > ptr && type2.length() > ptr) {

			if ((int)type1[ptr] < (int)type2[ptr]) {
				return;
			}
			else if ((int)type1[ptr] == (int)type2[ptr]) {
				// Do nothing, move ptr to the right
			}
			else {
				//printf("Swapping %d %d   ", (int) type1[ptr], (int) type2[ptr]);
				//cout << type1[ptr] << "   " << type2[ptr] << endl;
				//cout << type1 << "    " << type2 << endl << endl;
				swap(type1, type2);
				return;
			}
			ptr++;
		}
		if (type1.length() > type2.length()) {
			swap(type1, type2);
		}
		//printf("\n");
	}

	enum STATE { INACTIVE, BONDTYPES, ANGLETYPES, DIHEDRALTYPES, PAIRTYPES };
	static STATE setState(string s, STATE current_state) {
		if (s == "bondtypes" || s == "bonds")
			return BONDTYPES;
		if (s == "angletypes" || s == "angles")
			return ANGLETYPES;
		if (s == "dihedraltypes")
			return DIHEDRALTYPES;
		if (s == "pairs")
			return PAIRTYPES;
		return current_state;
	}
	static vector<Bondtype> parseFFBondtypes(vector<vector<string>> rows) {
		STATE current_state = INACTIVE;

		vector<Bondtype> records;

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
				records.push_back(Bondtype(row[0], row[1], stof(row[3]), stof(row[4])));
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

	static vector<Bondtype> parseTopolBondtypes(vector<vector<string>> rows) {
		STATE current_state = INACTIVE;
		vector<Bondtype> records;

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
				records.push_back(Bondtype(stoi(row[0]), stoi(row[1])));
				break;
			default:
				break;
			}

			if (current_state == PAIRTYPES)
				break;
		}
		printf("%d bonds found in topology file\n", records.size());
		return records;	
	}

	static void assignTypesFromAtomIDs(vector<Bondtype>* topol_bonds, vector<Atom> atoms) {
		for (int i = 0; i < topol_bonds->size(); i++) {
			Bondtype* bond = &topol_bonds->at(i);


			bond->type1 = atoms.at(bond->id1 - 1).atomtype_bond;	// Minus 1 becuase the bonds type1 is 1-indexed, and atoms vector is 0 indexed
			bond->type2 = atoms.at(bond->id2 - 1).atomtype_bond;
			bond->sort();
			//cout << bond->type1 << '\t' << bond->type2 << endl;;
		}
	}




	static Bondtype getBondFromTypes(Bondtype* query_type, vector<Bondtype>* FF_bondtypes) { 		
		float best_likeness = 0;
		Bondtype best_bond;
		for (Bondtype bond : *FF_bondtypes) {
			float likeness = calcLikeness(query_type->type1, bond.type1) * calcLikeness(query_type->type2, bond.type2);
			if (likeness > best_likeness) {
				best_likeness = likeness;
				best_bond = bond;
			}
		}
		if (best_likeness > 0.01f)
			return best_bond;

		cout << "Failed to match bond types.\n Closest match " << best_bond.type1 << "    " << best_bond.type2;
		printf("Likeness %f\n", best_likeness);		
		printf("Topol ids: %d %d\n", query_type->id1, query_type->id2);
		cout << query_type->type1 << '\t' << query_type->type2 << endl;
		exit(0);
	}

	static void assignFFParametersFromBondtypes(vector<Bondtype>* topol_bonds, vector<Bondtype>* FF_bondtypes) {
		for (int i = 0; i < topol_bonds->size(); i++) {
			if (!(i%10))
			printf("\rAssigning FF parameters to bond %06d of %06d", i, topol_bonds->size());
			Bondtype* bond = &topol_bonds->at(i);

			Bondtype appropriateForcefield = getBondFromTypes(bond, FF_bondtypes);	// This might not return the correct one, as it tries to fit the atomtypes_bond to atomtypes_bond known in the CHARMM forcefield
			
			bond->kb = appropriateForcefield.kb;
			bond->b0 = appropriateForcefield.b0;
		}
		printf("\n");
	}
};




struct Angletype {
	Angletype() {}
	Angletype(string t1, string t2, string t3) : type1(t1), type2(t2), type3(t3) {
		sort();
	}
	Angletype(string t1, string t2, string t3, float t0, float kt) : type1(t1), type2(t2), type3(t3), theta0(t0), ktheta(kt) {
		sort();
	}
	Angletype(int id1, int id2, int id3) : id1(id1), id2(id2), id3(id3) {
	}

	string type1, type2, type3;			// left, middle, right
	int id1, id2, id3;			// bonds from .top only has these values! 
	float theta0;
	float ktheta;


	void sort(string* leftmost, string* rightmost) {
		int ptr = 0;
		while (leftmost->length() > ptr && rightmost->length() > ptr) {
			if ((int)(*leftmost)[ptr] < (int)(*rightmost)[ptr]) {
				return;
			}
			else if ((int)(*leftmost)[ptr] > (int)(*rightmost)[ptr]) {
				swap(*leftmost, *rightmost);
				return;
			}
			ptr++;
		}
		if (leftmost->length() > rightmost->length()) {
			swap(*leftmost, *rightmost);
		}
	}
	void sort() {
		sort(&type1, &type3);		// We can ONLY sort these two, as type2 MUSTTT be in the middle!!
		//sort(&type1, &type2);
		//sort(&type2, &type3);
	}


	enum STATE { INACTIVE, BONDTYPES, ANGLETYPES, DIHEDRALTYPES };
	static STATE setState(string s, STATE current_state) {
		if (s == "bondtypes" || s == "bonds")
			return BONDTYPES;
		if (s == "angletypes" || s == "angles")
			return ANGLETYPES;
		if (s == "dihedraltypes" || s == "dihedrals")
			return DIHEDRALTYPES;
		return current_state;
	}

	static vector<Angletype> parseFFAngletypes(vector<vector<string>> rows) {
		STATE current_state = INACTIVE;
		vector<Angletype> angletypes;

		for (vector<string> row : rows) {

			if (row.size() == 3) {
				current_state = setState(row[1], current_state);
				continue;
			}



			switch (current_state)
			{
			case INACTIVE:
				break;
			case ANGLETYPES:
				if (row.size() < 8) {
					continue;													// TODO: dont just run from your problems, man...
				}
				angletypes.push_back(Angletype(row[0], row[1], row[2], stof(row[4]), stof(row[5])));
				break;
			case DIHEDRALTYPES:
				break;
			default:
				break;
			}


			if (current_state == DIHEDRALTYPES)
				break;
		}
		printf("%d angletypes in forcefield\n", angletypes.size());
		return angletypes;
	}

	static vector<Angletype> parseTopolAngletypes(vector<vector<string>> rows) {
		STATE current_state = INACTIVE;
		vector<Angletype> records;

		for (vector<string> row : rows) {
			if (row.size() == 3) {
				if (row[0][0] == '[') {
					current_state = setState(row[1], current_state);
					continue;
				}
			}



			switch (current_state)
			{
			case ANGLETYPES:
				records.push_back(Angletype(stoi(row[0]), stoi(row[1]), stoi(row[2])));
				break;
			default:
				break;
			}


			if (current_state == DIHEDRALTYPES)
				break;
		}
		printf("%d angles found in topology file\n", records.size());
		return records;
	}

	static void assignTypesFromAtomIDs(vector<Angletype>* topol_angles, vector<Atom> atoms) {
		for (int i = 0; i < topol_angles->size(); i++) {
			Angletype* angle = &topol_angles->at(i);

			angle->type1 = atoms.at(angle->id1 - 1).atomtype_bond;	// Minus 1 becuase the bonds type1 is 1-indexed, and atoms vector is 0 indexed
			angle->type2 = atoms.at(angle->id2 - 1).atomtype_bond;
			angle->type3 = atoms.at(angle->id3 - 1).atomtype_bond;
			angle->sort();
			//cout << bond->type1 << '\t' << bond->type2 << endl;;
		}
	}



	static Angletype getAngleFromTypes(Angletype* query_type, vector<Angletype>* FF_angletypes) {
		float best_likeness = 0;
		Angletype best_angle;
		for (Angletype angle : *FF_angletypes) {
			float likeness = calcLikeness(query_type->type1, angle.type1) * calcLikeness(query_type->type2, angle.type2) * calcLikeness(query_type->type3, angle.type3);
			if (likeness > best_likeness) {
				best_likeness = likeness;
				best_angle = angle;
			}
			if (likeness == 1)
				break;
		}
		if (best_likeness > 0.01f)
			return best_angle;

		cout << "\n\n\nFailed to match angle types.\n Closest match " << best_angle.type1 << "    " << best_angle.type2 << "    " << best_angle.type3 << endl;
		printf("Likeness %f\n", best_likeness);
		printf("Topol ids: %d %d %d\n", query_type->id1, query_type->id2, query_type->id3);
		cout << query_type->type1 << '\t' << query_type->type2 << '\t' << query_type->type3 << endl;
		exit(0);
	}
	
	static void assignFFParametersFromAngletypes(vector<Angletype>* topol_angles, vector<Angletype>* FF_angletypes) {
		for (int i = 0; i < topol_angles->size(); i++) {
			printf("\rAssigning FF parameters to angle %06d of %06d", i, topol_angles->size());

			Angletype* angle = &topol_angles->at(i);

			Angletype appropriateForcefield = getAngleFromTypes(angle, FF_angletypes);	// This might not return the correct one, as it tries to fit the atomtypes_bond to atomtypes_bond known in the CHARMM forcefield

			angle->theta0 = appropriateForcefield.theta0;
			angle->ktheta = appropriateForcefield.ktheta;
		}
		printf("\n");

	}
	
};




struct Dihedraltype {
	Dihedraltype() {}
	Dihedraltype(string t1, string t2, string t3, string t4) : type1(t1), type2(t2), type3(t3), type4(t4) {
		sort();
	}
	Dihedraltype(string t1, string t2, string t3, string t4, float phi0, float kphi) : type1(t1), type2(t2), type3(t3), type4(t4), phi0(phi0), kphi(kphi) {
		sort();
	}
	Dihedraltype(int id1, int id2, int id3, int id4) : id1(id1), id2(id2), id3(id3), id4(id4) {
	}

	string type1, type2, type3, type4;			// left, lm, rm, right
	int id1, id2, id3, id4;			// bonds from .top only has these values! 
	float phi0;
	float kphi;


	void flip() {
		swap(type1, type4);
		swap(type2, type3);
	}
	void sort() {		
		if (type1 != type4) {
			if (!isSorted(&type1, &type4)) {
				flip();
			}
		}
		else {			// In case the outer two is identical, we check the inner two.
			if (!isSorted(&type2, &type3)) {
				flip();
			}
		}
	}


	enum STATE { INACTIVE, BONDTYPES, ANGLETYPES, DIHEDRALTYPES, CMAP };
	static STATE setState(string s, STATE current_state) {
		if (s == "bondtypes" || s == "bonds")
			return BONDTYPES;
		if (s == "angletypes" || s == "angles")
			return ANGLETYPES;
		if (s == "dihedraltypes" || s == "dihedrals")
			return DIHEDRALTYPES;
		if (s == "cmap")
			return CMAP;
		return current_state;
	}

	static vector<Dihedraltype> parseFFDihedraltypes(vector<vector<string>> rows) {
		STATE current_state = INACTIVE;
		vector<Dihedraltype> dihedraltypes;

		for (vector<string> row : rows) {

			if (row.size() == 3) {
				current_state = setState(row[1], current_state);
				continue;
			}



			switch (current_state)
			{
			case DIHEDRALTYPES:
				if (row.size() < 7) {
					continue;													// TODO: dont just run from your problems, man...
				}
				dihedraltypes.push_back(Dihedraltype(row[0], row[1], row[2], row[3], stof(row[5]), stof(row[6])));
				break;
			default:
				break;
			}

		}
		printf("%d dihedraltypes in forcefield\n", dihedraltypes.size());
		return dihedraltypes;
	}

	static vector<Dihedraltype> parseTopolDihedraltypes(vector<vector<string>> rows) {
		STATE current_state = INACTIVE;
		vector<Dihedraltype> records;

		for (vector<string> row : rows) {
			if (row.size() == 3) {
				if (row[0][0] == '[') {
					current_state = setState(row[1], current_state);
					continue;
				}
			}



			switch (current_state)
			{
			case DIHEDRALTYPES:
				records.push_back(Dihedraltype(stoi(row[0]), stoi(row[1]), stoi(row[2]), stoi(row[3])));
				break;
			default:
				break;
			}
		}
		printf("%d dihedrals found in topology file\n", records.size());
		return records;
	}

	static void assignTypesFromAtomIDs(vector<Dihedraltype>* topol_dihedrals, vector<Atom> atoms) {
		//printf("GOt here!\n");
		for (int i = 0; i < topol_dihedrals->size(); i++) {
			Dihedraltype* dihedral = &topol_dihedrals->at(i);

			//printf("Accessing atoms %d %d %d %d\n", dihedral->id1, dihedral->id2, dihedral->id3, dihedral->id4);
			dihedral->type1 = atoms.at(dihedral->id1 - 1).atomtype_bond;	// Minus 1 becuase the bonds type1 is 1-indexed, and atoms vector is 0 indexed
			dihedral->type2 = atoms.at(dihedral->id2 - 1).atomtype_bond;
			dihedral->type3 = atoms.at(dihedral->id3 - 1).atomtype_bond;
			dihedral->type4 = atoms.at(dihedral->id4 - 1).atomtype_bond;
			dihedral->sort();
			//cout << bond->type1 << '\t' << bond->type2 << endl;;
		}
	}



	static Dihedraltype getDihedralFromTypes(Dihedraltype* query_type, vector<Dihedraltype>* forcefield) {
		float best_likeness = 0;
		Dihedraltype best_match;
		for (Dihedraltype dihedral : *forcefield) {
			float likeness = calcLikeness(query_type->type1, dihedral.type1) * calcLikeness(query_type->type2, dihedral.type2) * calcLikeness(query_type->type3, dihedral.type3) * calcLikeness(query_type->type4, dihedral.type4);
			if (likeness > best_likeness) {
				best_likeness = likeness;
				best_match = dihedral;
			}
			if (likeness == 1)
				break;
		}
		if (best_likeness > 0.001f)
			return best_match;

		cout << "\n\n\nFailed to match angle types.\nClosest match: " << best_match.type1 << "    " << best_match.type2 << "    " << best_match.type3 << "    " << best_match.type4 << endl;
		printf("Likeness %f\n", best_likeness);
		printf("Topol ids: %d %d %d %d\n", query_type->id1, query_type->id2, query_type->id3, query_type->id4);
		cout << query_type->type1 << '\t' << query_type->type2 << '\t' << query_type->type3 << '\t' << query_type->type4 << endl;
		exit(0);
	}

	static void assignFFParametersFromDihedraltypes(vector<Dihedraltype>* topol_dihedrals, vector<Dihedraltype>* forcefield) {
		for (int i = 0; i < topol_dihedrals->size(); i++) {
			printf("\rAssigning FF parameters to dihedral %06d of %06d", i, topol_dihedrals->size());

			Dihedraltype* dihedral = &topol_dihedrals->at(i);

			Dihedraltype appropriateForcefield = getDihedralFromTypes(dihedral, forcefield);	// This might not return the correct one, as it tries to fit the atomtypes_bond to atomtypes_bond known in the CHARMM forcefield

			dihedral->phi0 = appropriateForcefield.phi0;
			dihedral->kphi = appropriateForcefield.kphi;
		}
		printf("\n");
	}
};








