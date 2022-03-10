#include "CompoundBuilder.h"



float massFromAtom(char atom) {
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
	case 'H':
		return 1.f;
	default:
		printf("Element %c not found", atom);
		exit(1);
	}
}



Compound CompoundBuilder::buildMolecule(string pdb_path, string itp_path, int max_residue_id)
{
	//vector<vector<string>> pdb_data = readFile(pdb_path);
	//vector<Record_ATOM> atom_data = parsePDB(pdb_path);
	vector<Record_ATOM> atom_data = parseGRO(pdb_path);
	vector<vector<string>> top_data = parseTOP(itp_path);


	Compound compound;													// position of particles stored in compund.particles.pos_tsub1

	loadParticles(&compound, &atom_data, max_residue_id, true);
	printf("%d particles added\n", compound.n_particles);
	loadTopology(&compound, &top_data, particle_id_map);
	printf("%d pairbonds added\n", compound.n_pairbonds);
	printf("%d anglebonds added\n", compound.n_anglebonds);

	calcParticleSphere(&compound);



	delete[] particle_id_map;
	return compound;
}

void CompoundBuilder::loadParticles(Compound* compound, vector<CompoundBuilder::Record_ATOM>* pdb_data, int max_residue_id, bool ignore_protons)
{
	int max_atom_cnt = 1 << 24;
	particle_id_map = new int[max_atom_cnt];
	for (int i = 0; i < max_atom_cnt; i++)
		particle_id_map[i] = -1;


	for (Record_ATOM record : *pdb_data) {

			/*
			int particle_id = stoi(record[1]);
			string atom_name = record[2];
			string monomer_name = record[3];
			cout << record[3] <<endl;
			cout << record[4] << endl;
			int monomer_id = stoi(record[4]);
			Float3 coord(stod(record[5]), stod(record[6]), stod(record[7]));
			coord *= 0.1f;	// convert A to nm
			*/
			//printf("res %d   max res %d\n", record.residue_seq_number, max_residue_id);
			if (record.residue_seq_number > max_residue_id)
				break;

			if (ignore_protons && record.atom_name[0] == 'H')
				continue;


			particle_id_map[record.atom_serial_number] = compound->n_particles;


			//if (stoi(record[4]) > max_monomer_id)
				//break;


			//particle_id_map[particle_id] = compound->n_particles;
			float mass = massFromAtom(record.atom_name[0]) / 1000.f;							// kg/mol;

			compound->atom_types[compound->n_particles] = FFM.atomTypeToIndex(record.atom_name[0]);
			compound->particles[compound->n_particles] = CompactParticle(record.coordinate);
			compound->n_particles++;
	}
}

void CompoundBuilder::loadTopology(Compound* compound, vector<vector<string>>* top_data, int* particle_id_map)
{
	TopologyMode mode = INACTIVE;
	for (vector<string> record : *top_data) {
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

void CompoundBuilder::calcParticleSphere(Compound* compound) {
	Float3 com = calcCOM(compound);

	float furthest = LONG_MIN;
	float closest = LONG_MAX;
	int closest_index = 0;

	for (int i = 0; i < compound->n_particles; i++) {
		float dist = (compound->particles[i].pos_tsub1 - com).len();

		if (dist < closest) {
			closest = dist;
			closest_index = i;
		}
		furthest = max(furthest, dist);
	}

	compound->key_particle_index = closest_index;
	compound->confining_particle_sphere = furthest;
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

Float3 CompoundBuilder::calcCOM(Compound* compound) {
	Float3 com(0.f);
	for (int i = 0; i < compound->n_particles; i++) {
		com += compound->particles[i].pos_tsub1;
	}
	return com * (1.f / (float) compound->n_particles);
}

vector<vector<string>> CompoundBuilder::parseTOP(string path)		// Naive file segmentation, DANGEROUS!
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

vector<CompoundBuilder::Record_ATOM> CompoundBuilder::parsePDB(string path)
{
	fstream file;
	file.open(path);
	int line_cnt = 0;

	int endpoints[] = { 4, 11, 16 , 17, 20, 22, 26, 27, 38, 46, 54 };


	vector<Record_ATOM> records;

	string line;
	while (getline(file, line)) {
		stringstream ss(line);
		string row_type;
		getline(ss, row_type, ' ');
		if (row_type != "ATOM")
			continue;

		vector<string> data_buffer;


		int ptr = 0;
		for (int stop : endpoints) {
			string word = "";
			
			while (ptr < stop) {
				if (line[ptr] != ' ')
					word = word + line[ptr];
				ptr++;
			}
			//cout << "Word:" << word << endl;
			data_buffer.push_back(word);
		}

		cout << data_buffer[6] << endl;
		//printf("%d\n", stoi(data_buffer[6]));
		
		records.push_back(Record_ATOM(
			stoi(data_buffer[1]),
			data_buffer[2],
			data_buffer[3][0],			// alt loc
			data_buffer[4],
			data_buffer[5][0],		// chain id
			stoi(data_buffer[6]),
			data_buffer[7][0],
			Float3(stod(data_buffer[8]), stod(data_buffer[9]), stod(data_buffer[10])) * 0.1	// Convert A to nm right off the bat!
		));
		
	}
	return records;
}

vector<CompoundBuilder::Record_ATOM> CompoundBuilder::parseGRO(string path)
{
	fstream file;
	file.open(path);
	int line_cnt = 0;
	string line;
	string space_delimiter = " ";

	vector<Record_ATOM> records;

	while (getline(file, line)) {
		vector<string> words;
		stringstream ss(line);
		string word;
		while (getline(ss, word, ' ')) {
			if (word == "" || word == "[" || word == "]") {
				continue;
			}
			words.push_back(word);
		}


		ResidueComboId res = parseResidueID(words.at(0));

		if (!res.valid || words.size() != 6)
			continue;

		Record_ATOM record(
			stoi(words.at(2)),
			words.at(1),
			' ',
			res.name,
			' ',
			res.id,
			' ',
			Float3(stod(words.at(3)), stod(words.at(4)), stod(words.at(5)))
		);
		records.push_back(record);
	}

	

	return records;
}

CompoundBuilder::ResidueComboId CompoundBuilder::parseResidueID(string s)
{
	string id = "";
	string name = "";
	for (char c : s) {
		if (isAsciiNumber(c))
			id = id + c;
		else
			name = name + c;
	}
	
	if (id == "" || name == "")
		return ResidueComboId();
	return ResidueComboId(stoi(id), name);
}

bool CompoundBuilder::isAsciiNumber(char c)
{
	return (c > 47 && c < 58);
}

