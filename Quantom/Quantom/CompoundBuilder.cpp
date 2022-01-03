#include "CompoundBuilder.h"



float mass_from_atom(char atom) {
	switch (atom)
	{
	case 'C':
		return 12.0096; //g/mol
	case 'O':
		return 15.9994;	
	case 'P':
		return 30.9738;
	case 'N':
		return 14.0067;
	default:
		printf("Element %c not found", atom);
		exit(1);
	}
}



Compound CompoundBuilder::buildMolecule(string pdb_path, string itp_path)
{
	vector<vector<string>> pdb_data = readFile(pdb_path);
	vector<vector<string>> itp_data = readFile(itp_path);


	Compound compound;

	loadParticles(&compound, &pdb_data, 1);
	printf("%d particles added\n", compound.n_particles);
	loadTopology(&compound, &itp_data, particle_id_map);
	printf("%d pairbonds added\n", compound.n_pairbonds);
	printf("%d anglebonds added\n", compound.n_anglebonds);



	delete[] particle_id_map;
	return compound;
}

void CompoundBuilder::loadParticles(Compound* compound, vector<vector<string>>* pdb_data, int max_monomer_id)
{
	int size = 1 << 24;
	particle_id_map = new int[size];
	for (int i = 0; i < size; i++)
		particle_id_map[i] = -1;


	for (vector<string> record : *pdb_data) {
		if (record[0] == "ATOM") {

			
			int particle_id = stoi(record[1]);
			string atom_name = record[2];
			string monomer_name = record[3];
			int monomer_id = stoi(record[4]);
			Float3 coord(stod(record[5]), stod(record[6]), stod(record[7]));
			coord *= 0.1f;	// convert A to nm

			if (stoi(record[4]) > max_monomer_id)
				break;


			particle_id_map[particle_id] = compound->n_particles;
			float mass = mass_from_atom(atom_name[0]) / 1000.f;							// kg/mol;


			compound->particles[compound->n_particles++] = CompactParticle(mass, coord);









			//if (compound->n_particles == 5)
				//break;
		}
	}
}

void CompoundBuilder::loadTopology(Compound* compound, vector<vector<string>>* itp_data, int* particle_id_map)
{
	TopologyMode mode = INACTIVE;
	for (vector<string> record : *itp_data) {
		if (record.size() == 0) {
			mode = INACTIVE;
			continue;
		}
		
		if (mode == INACTIVE) {
			mode = setMode(record[0]);
			continue;
		}

		addGeneric(compound, &record, mode);
	}
}

CompoundBuilder::TopologyMode CompoundBuilder::setMode(string entry)
{
	if (entry == "bonds")
		return BOND;
	if (entry == "angles")
		return ANGLE;
	if (entry == "dihedrals")
		return DIHEDRAL;
	return INACTIVE;
}

void CompoundBuilder::addGeneric(Compound* compound, vector<string>* record, TopologyMode mode)
{
	switch (mode)
	{
	case CompoundBuilder::INACTIVE:
		return;
	case CompoundBuilder::BOND:
		return addBond(compound, record);
	case CompoundBuilder::ANGLE:
		return addAngle(compound, record);
	case CompoundBuilder::DIHEDRAL:
		return addDihedral(compound, record);
	default:
		return;
	}
}

void CompoundBuilder::addBond(Compound* compound, vector<string>* record)
{
	int particle_indexes[2] = { particle_id_map[stoi((*record)[0])],
							particle_id_map[stoi((*record)[1])],
	};	// left, right

	for (int i = 0; i < 2; i++) {
		if (particle_indexes[i] == -1) // In this case, either of the atoms have not been loaded, thus we cannot bind them
			return;
	}

	double dist = (compound->particles[particle_indexes[0]].pos_tsub1 - compound->particles[particle_indexes[1]].pos_tsub1).len();
	//printf("Calculated ref dist: %f\n", dist);
	compound->pairbonds[compound->n_pairbonds++] = PairBond(dist, particle_indexes[0], particle_indexes[1]);
	//printf("Pairbond error: %f\n", OH_refdist - ();
}

void CompoundBuilder::addAngle(Compound* compound, vector<string>* record)
{
	int particle_indexes[3] = {	particle_id_map[stoi((*record)[0])],
						particle_id_map[stoi((*record)[1])],
						particle_id_map[stoi((*record)[2])]
	};	// left, middle, right
	
	for (int i = 0; i < 3; i++) {
		if (particle_indexes[i] == -1) 
			return;
	}

	double angle = Float3::getAngle(compound->particles[particle_indexes[0]].pos_tsub1, compound->particles[particle_indexes[1]].pos_tsub1, compound->particles[particle_indexes[2]].pos_tsub1);
	compound->anglebonds[compound->n_anglebonds++] = AngleBond(angle, particle_indexes[0], particle_indexes[1], particle_indexes[2]);
}

void CompoundBuilder::addDihedral(Compound* compound, vector<string>* record)
{
}

vector<vector<string>> CompoundBuilder::readFile(string path)
{
	fstream file;
	file.open(path);
	int line_cnt = 0;
	string line;
	string space_delimiter = " ";

	vector<vector<string>> records;

	while (getline(file, line)) {
		vector<string> record;


		stringstream ss(line);
		string word;
		while (getline(ss, word, ' ')) {
			if (word == "" || word == "[" || word == "]") {
				continue;
			}
			record.push_back(word);
		}
		if (record.size() > 0) {
			if (record[0] == ";")
				continue;
		}
		records.push_back(record);
	}

	for (auto record : records) {
		for (auto w : record) {
			//cout << w << " ";
		}
		//printf("\n");
	}
	return records;
}
