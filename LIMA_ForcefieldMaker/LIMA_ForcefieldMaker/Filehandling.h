#pragma once

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>

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





