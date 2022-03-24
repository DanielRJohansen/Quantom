#pragma once


#include "Forcefield.cuh"





void ForceFieldMaker::buildForcefield() {
	vector<vector<string>> summary_rows = Filehandler::readFile("C:\\Users\\Daniel\\git_repo\\Quantom\\ForcefieldSummary.txt", INT_MAX, true);
	vector<vector<string>> forcefield_rows = Filehandler::readFile("C:\\Users\\Daniel\\git_repo\\Quantom\\Forcefield.txt", INT_MAX, true);

	NBAtomtype* nb_atomtypes = parseAtomTypes(summary_rows);					// 1 entry per type in compressed forcefield
	int* nb_atomtype_ids = parseAtomTypeIDs(forcefield_rows);				// 1 entry per atom in conf

	PairBond* bonds = parseBonds(forcefield_rows);
	AngleBond* angles = parseAngles(forcefield_rows);
	DihedralBond* dihedrals = parseDihedrals(forcefield_rows);

	//exit(0);
}




NBAtomtype* ForceFieldMaker::parseAtomTypes(vector<vector<string>> summary_rows) {
	NBAtomtype* atomtypes = new NBAtomtype[10000];
	int ptr = 0;
	STATE current_state = INACTIVE;

	for (vector<string> row : summary_rows) {
		if (newParseTitle(row)) {
			current_state = setState(row[1], current_state);
			continue;
		}

		

		if (current_state == FF_NONBONDED) {
			for (string e : row)
				cout << e << '\t';
			printf("\n");
			atomtypes[ptr++] = NBAtomtype(stof(row[2]), stof(row[3]), stof(row[4]));
		}			
	}
	printf("%d NB_Atomtypes loaded\n", ptr);
}

int* ForceFieldMaker::parseAtomTypeIDs(vector<vector<string>> forcefield_rows) {	// returns the nonbonded atomtype
	int* atomtype_ids = new int[10000];
	int ptr = 0;
	STATE current_state = INACTIVE;

	for (vector<string> row : forcefield_rows) {
		if (newParseTitle(row)) {
			current_state = setState(row[1], current_state);
			continue;
		}

		if (current_state == NB_ATOMTYPES)
			atomtype_ids[ptr++] = stoi(row[1]);
	}
	printf("%d NB_Atomtype_IDs loaded\n", ptr);
}

PairBond* ForceFieldMaker::parseBonds(vector<vector<string>> forcefield_rows) {
	PairBond* bonds = new PairBond[10000];
	int ptr = 0;
	STATE current_state = INACTIVE;

	for (vector<string> row : forcefield_rows) {
		if (newParseTitle(row)) {
			current_state = setState(row[1], current_state);
			continue;
		}

		if (current_state == BONDS) {
			bonds[ptr++] = PairBond(stoi(row[0]), stoi(row[1]), stof(row[4]), stof(row[5]));
		}
			
	}
	printf("%d bonds loaded\n", ptr);
	return bonds;
}

AngleBond* ForceFieldMaker::parseAngles(vector<vector<string>> forcefield_rows) {
	AngleBond* angles = new AngleBond[10000];
	int ptr = 0;
	STATE current_state = INACTIVE;

	for (vector<string> row : forcefield_rows) {
		if (newParseTitle(row)) {
			current_state = setState(row[1], current_state);
			continue;
		}

		if (current_state == ANGLES) {
			angles[ptr++] = AngleBond(stof(row[6]), stof(row[7]));
		}

	}
	printf("%d angles loaded\n", ptr);
	return angles;
}

DihedralBond* ForceFieldMaker::parseDihedrals(vector<vector<string>> forcefield_rows) {
	DihedralBond* dihedrals = new DihedralBond[10000];
	int ptr = 0;
	STATE current_state = INACTIVE;

	for (vector<string> row : forcefield_rows) {
		if (newParseTitle(row)) {
			current_state = setState(row[1], current_state);
			continue;
		}

		if (current_state == DIHEDRALS) {
			dihedrals[ptr++] = DihedralBond(stof(row[8]), stof(row[9]));
		}

	}
	printf("%d dihdrals loaded\n", ptr);
	return dihedrals;
}
