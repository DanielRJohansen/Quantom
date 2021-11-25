#include "CompoundBuilder.h"







Molecule1 CompoundBuilder::buildMolecule(string pdb_path)
{
	int n_monomers = countMonomers(pdb_path);



	


	string line;
	//getline(file, line);


	
	

	exit(0);


	return Molecule1();
}

int CompoundBuilder::countMonomers(string path)
{
	fstream file;
	file.open(path);
	string line;

	while (getline(file, line)) {
		cout << line << endl;
		for (int i = 0; i < line.size(); i++)
			cout << line[i] << endl;

		exit(1);
	}
	return 0;
}
