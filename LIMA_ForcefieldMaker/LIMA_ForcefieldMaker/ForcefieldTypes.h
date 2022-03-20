#pragma once

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>



using namespace std;



struct FF_nonbonded {
	FF_nonbonded() {}
	FF_nonbonded(string t) : type(t){}		// This is for loading form the conf file, for easy comparisons
	FF_nonbonded(string t, int a, double m, double s, double e) : type(t), atnum(a), mass(m), sigma(s), epsilon(e) {}
	// Official parameters
	string type;
	int atnum;
	double mass, sigma, epsilon;

	// LIMA parameters
	bool is_present_in_simulation = false;
	int simulation_specific_id;



	bool operator==(const FF_nonbonded a) { return (type == a.type); }







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
					//printf("size %d\n", row.size());
					//for (int i = 0; i < row.size(); i++)
					//	cout << "HERE: " << row[i] << endl;
					//current_state = INACTIVE;
					//goto START;
				}
				//cout << row[0] << endl;
				records.push_back(FF_nonbonded(row[0], stoi(row[1]), stod(row[2]), stod(row[5]), stod(row[6])));
				break;
			case PAIRTYPES:
				break;
			default:
				break;
			}
		}
		return records;
	}
	static vector<FF_nonbonded> parseConf(vector<vector<string>> rows) {		// KNOWN ERROR HERE. SOMETIMES SPACES IN FILES ARE TABS, AND IT DOESN'T RECOGNISE A ROW AS A PROPER ENTRY!!
		STATE current_state = INACTIVE;

		vector<FF_nonbonded> records;

		for (vector<string> row : rows) {
			if (row.size() != 6) {
				continue;
			}

			bool type_already_found = false;
			for (FF_nonbonded record : records) {
				if (record.type == row[1]) {
					type_already_found = true;
					break;
				}
			}
			if (!type_already_found)
				records.push_back(FF_nonbonded(row[1]));
		}

		printf("%d atom types in conf\n", records.size());

		return records;
	}
};




