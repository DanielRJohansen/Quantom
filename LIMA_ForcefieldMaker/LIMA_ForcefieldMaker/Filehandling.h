#pragma once

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>

#include "ForcefieldTypes.h"

using namespace std;


bool ignoreRow(vector<char> ignores, string line) {
	if (line.length() == 0)
		return true;
	//if (*find(ignores.begin(), ignores.end(), line[0]) == line[0])
		//return true;
	for (char c : ignores) {
		if (line[0] == c)
			return true;
	}
	return false;
}

vector<vector<string>> readFile(string path, int end_at = INT_MAX) {
	fstream file;
	file.open(path);


	vector<char> ignores = { ';', '#' };





	vector<vector<string>> rows;
	int row_cnt = 0;
	int ignore_cnt = 0;

	string line;
	while (getline(file, line)) {



		if (ignoreRow(ignores, line)) {
			ignore_cnt++;
			continue;
		}

		vector<string> row;
		stringstream ss(line);
		string element;
		while (getline(ss, element, ' ')) {
			if (element.length() > 0)
				row.push_back(element);
		}


		rows.push_back(row);
		row_cnt++;

		if (row_cnt >= end_at)
			break;
	}
	printf("%d rows read. %d rows ignored\n", row_cnt, ignore_cnt);

	return rows;
}





struct FFOutHelpers{
	string includes = "#pragma once\n";

	static string titleH1(string text) {
		return "// ----------------################ " + text + " ################---------------- //\n\n";
	}
	static string titleH2(string text) {
		return "// ---------------- " + text + " ---------------- //\n";
	}
	static string titleH3(string text) {
		return "// " + text + " //\n";
	}
	static string parserTitle(string text) {
		return "# " + text + '\n';
	}
	static string endBlock() {
		return "\n\n\n\n";
	}
};








void printForcefieldSummary(string path, vector<FF_nonbonded> records_nonbonded, Map* map) {
	ofstream file(path, ofstream::out);
	if (!file.is_open()) {
		printf(("Failed to open file\n"));
		exit(0);
	}

	file << FFOutHelpers().includes;
	file << FFOutHelpers::titleH1("Forcefield Non-bonded");
	file << FFOutHelpers::titleH2("Mappings");
	file << FFOutHelpers::parserTitle("mappings");
	for (int i = 0; i < map->n_mappings; i++) {
		file << map->mappings[i].left << ";" << map->mappings[i].right << endl;
	}
	file << FFOutHelpers::endBlock();

	file << FFOutHelpers::titleH2("Non-bonded parameters");
	file << FFOutHelpers::titleH2("{atom_type\tatnum\tmass [g/mol]\tsigma []\tepsilon []}");
	file << FFOutHelpers::parserTitle("ff_nonbonded");
	for (FF_nonbonded record : records_nonbonded) {
		file << record.type << ';' << to_string(record.atnum_local) << ';' << to_string(record.mass) << ';' << to_string(record.sigma) << ';' << to_string(record.epsilon) << endl;
	}
	file << FFOutHelpers::endBlock();


	file.close();
}

void printForcefield(string path, vector<Atom> atoms, vector<Bondtype> bonds, vector<Angletype> angles, vector<Dihedraltype> dihedrals) {

	ofstream file(path, ofstream::out);
	if (!file.is_open()) {
		printf(("Failed to open file\n"));
		exit(0);
	}

	file << FFOutHelpers::titleH1("Forcefield Non-bonded");
	file << FFOutHelpers::titleH2("Atoms {particle id [simulation specific]\tatomtype_id [simulation specific]}");
	file << FFOutHelpers::parserTitle("atoms");
	for (Atom atom : atoms) {
		file << to_string(atom.id) << ";" << to_string(atom.atomtype_id) << endl;
	}
	file << FFOutHelpers::endBlock();






	file << FFOutHelpers::titleH2("Bonds");
	file << FFOutHelpers::titleH3("{ID_p1 \t ID_p2 \t b_0 \t k_b}");
	file << FFOutHelpers::parserTitle("bonds");
	for (Bondtype bond : bonds) {
		file << to_string(bond.id1) << ';' << to_string(bond.id2) << ';'
			<< bond.type1 << ';' << bond.type2 << ';'
			<< to_string(bond.b0) << ';' << to_string(bond.kb) << endl;
	}
	file << FFOutHelpers::endBlock();



	file << FFOutHelpers::titleH2("Angles");
	file << FFOutHelpers::titleH3("{ID_p1 \t ID_p2 \t ID_p3 \t theta_0 \t k_theta}");
	file << FFOutHelpers::parserTitle("angles");
	for (Angletype angle : angles) {
		file << to_string(angle.id1) << ';' << to_string(angle.id2) << ';' << to_string(angle.id3) << ';'
			<< angle.type1 << ';' << angle.type2 << ';' << angle.type3 << ';'
			<< to_string(angle.theta0) << ';' << to_string(angle.ktheta) << endl;
	}
	file << FFOutHelpers::endBlock();



	file << FFOutHelpers::titleH2("Dihedrals");
	file << FFOutHelpers::titleH3("{ID_p1 \t ID_p2 \t ID_p3 \t ID_p4 \t phi_0 \t k_phi}");
	file << FFOutHelpers::parserTitle("dihedrals");
	for (Dihedraltype dihedral: dihedrals) {
		file << to_string(dihedral.id1) << ';' << to_string(dihedral.id2) << ';' << to_string(dihedral.id3) << ';' << to_string(dihedral.id4) << ';'
			<< dihedral.type1 << ';' << dihedral.type2 << ';' << dihedral.type3 << ';' << dihedral.type4 << ';'
			<< to_string(dihedral.phi0) << ';' << to_string(dihedral.kphi) << endl;
	}
	file << FFOutHelpers::endBlock();

	file.close();
}

