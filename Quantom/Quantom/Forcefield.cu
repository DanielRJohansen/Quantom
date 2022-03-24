#pragma once


#include "Forcefield.cuh"





void ForceFieldMaker::buildForcefield() {
	vector<vector<string>> summary_rows = Filehandler::readFile("C:\\Users\\Daniel\\git_repo\\Quantom\\ForcefieldSummary.txt", INT_MAX, true);
	vector<vector<string>> forcefield_rows = Filehandler::readFile("C:\\Users\\Daniel\\git_repo\\Quantom\\Forcefield.txt", INT_MAX, true);

	NBAtomtype* nb_atomtypes = parseAtomTypes(summary_rows);					// 1 entry per type in compressed forcefield
	int* nb_atomtype_ids = parseAtomTypeIDs(forcefield_rows);				// 1 entry per atom in conf

	topol_bonds = parseBonds(forcefield_rows);
	topol_angles = parseAngles(forcefield_rows);
	topol_dihedrals = parseDihedrals(forcefield_rows);

	//exit(0);
}



PairBond* ForceFieldMaker::getBondType(int id1, int id2) {
	for (int i = 0; i < n_topol_bonds; i++) {
		if (topol_bonds[i].atom_indexes[0] == id1 && topol_bonds[i].atom_indexes[1] == id2) {
			return &topol_bonds[i];
		}
	}
	printf("Bond not found with ids %d %d\n", id1, id2);
	exit(0);
}

AngleBond* ForceFieldMaker::getAngleType(int id1, int id2, int id3) {
	for (int i = 0; i < n_topol_angles; i++) {
		if (topol_angles[i].atom_indexes[0] == id1 && topol_angles[i].atom_indexes[1] == id2 && topol_angles[i].atom_indexes[2] == id3) {
			return &topol_angles[i];
		}
	}
	printf("Angle not found with ids %d %d %d\n", id1, id2, id3);
	exit(0);
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
	n_topol_bonds = ptr;
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
			angles[ptr++] = AngleBond(stoi(row[0]), stoi(row[1]), stoi(row[2]), stof(row[6]) * (2*PI)/360.f , stof(row[7]));		// Convert degree to radians here!
		}

	}
	n_topol_angles = ptr;
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
	n_topol_dihedrals = ptr;
	printf("%d dihdrals loaded\n", ptr);
	return dihedrals;
}
