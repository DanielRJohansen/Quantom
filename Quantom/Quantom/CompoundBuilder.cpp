#include "CompoundBuilder.h"





Molecule CompoundBuilder::buildMolecule(string pdb_path, string itp_path, int max_residue_id)
{
	//vector<vector<string>> pdb_data = readFile(pdb_path);
	//vector<Record_ATOM> atom_data = parsePDB(pdb_path);
	particle_id_maps = new IDMap[MAX_ATOMS];

	vector<Record_ATOM> atom_data = parseGRO(pdb_path);
	vector<vector<string>> top_data = parseTOP(itp_path);


	//Compound compound;													// position of particles stored in compund.particles.pos_tsub1
	Molecule molecule;

	//loadParticles(&compound, &atom_data, max_residue_id, true);
	loadParticles(&molecule, &atom_data, max_residue_id, true);
	//printf("%d particles added\n", compound.n_particles);
	loadTopology(&molecule, &top_data);



	countElements(&molecule);
	printf("Molecule built\n\n\n");

	for (int i = 0; i < molecule.n_compounds; i++) {
		molecule.compounds[i].calcParticleSphere();
	}


	delete[] particle_id_maps;
	return molecule;
}





void CompoundBuilder::loadParticles(Molecule* molecule, vector<CompoundBuilder::Record_ATOM>* pdb_data, int max_residue_id, bool ignore_protons) {


	int current_res_id = 0;
	int current_compound_id = 0;
	Compound* current_compound = &molecule->compounds[current_compound_id];
	int prev_atom_id;

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


			//particle_id_map[record.atom_serial_number] = compound->n_particles;
			if (record.residue_seq_number != current_res_id) {
				if (!current_compound->hasRoomForRes()) {

					molecule->compound_bridge.addParticle(current_compound->atom_types[current_compound->n_particles - 1], current_compound->particles[current_compound->n_particles - 1]);		// Add prev atom to the bridge
					molecule->compound_bridge.addParticle(FFM.atomTypeToIndex(record.atom_name[0]), CompactParticle(record.coordinate));														// And new atom to the bridge
					current_compound_id++;
					current_compound = &molecule->compounds[current_compound_id];
					molecule->n_compounds++;
				}
				current_res_id = record.residue_seq_number;
			}

			particle_id_maps[record.atom_serial_number] = IDMap(record.atom_serial_number, current_compound_id, molecule->compounds[current_compound_id].n_particles);
			current_compound->addParticle(FFM.atomTypeToIndex(record.atom_name[0]), CompactParticle(record.coordinate));
			//compound->addParticle(FFM.atomTypeToIndex(record.atom_name[0]), CompactParticle(record.coordinate));
	}
}

void CompoundBuilder::loadTopology(Molecule* molecule, vector<vector<string>>* top_data)
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

		addGeneric(molecule, &record, mode);
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

void CompoundBuilder::addGeneric(Molecule* molecule, vector<string>* record, TopologyMode mode) {
	switch (mode)
	{
	case CompoundBuilder::INACTIVE:
		return;
	case CompoundBuilder::BOND:
		return addBond(molecule, record);
	case CompoundBuilder::ANGLE:
		return addAngle(molecule, record);
	case CompoundBuilder::DIHEDRAL:
		return addDihedral(molecule, record);
	default:
		return;
	}
}

void CompoundBuilder::addBond(Molecule* molecule, vector<string>* record) {
	Compound* compound = &molecule->compounds[particle_id_maps[stoi((*record)[0])].compound_id];
	int particle_indexes[2] = { 
		particle_id_maps[stoi((*record)[0])].local_id,
		particle_id_maps[stoi((*record)[1])].local_id
	};	// left, right

	for (int i = 0; i < 2; i++) {
		if (particle_indexes[i] == -1) // In this case, either of the atoms have not been loaded, thus we cannot bind them
			return;
	}

	double dist = (compound->particles[particle_indexes[0]].pos_tsub1 - compound->particles[particle_indexes[1]].pos_tsub1).len();			// TODO: Remove, information is in forcefield!
	//printf("Calculated ref dist: %f\n", dist);
	compound->pairbonds[compound->n_pairbonds++] = PairBond(dist, particle_indexes[0], particle_indexes[1]);
	//printf("Pairbond error: %f\n", OH_refdist - ();
}

void CompoundBuilder::addAngle(Molecule* molecule, vector<string>* record) {
	Compound* compound = &molecule->compounds[particle_id_maps[stoi((*record)[0])].compound_id];
	int particle_indexes[3] = {	
		particle_id_maps[stoi((*record)[0])].local_id,
		particle_id_maps[stoi((*record)[1])].local_id,
		particle_id_maps[stoi((*record)[2])].local_id
	};	// left, middle, right
	
	for (int i = 0; i < 3; i++) {
		if (particle_indexes[i] == -1) 
			return;
	}

	double angle = Float3::getAngle(compound->particles[particle_indexes[0]].pos_tsub1, compound->particles[particle_indexes[1]].pos_tsub1, compound->particles[particle_indexes[2]].pos_tsub1);
	compound->anglebonds[compound->n_anglebonds++] = AngleBond(angle, particle_indexes[0], particle_indexes[1], particle_indexes[2]);
}

void CompoundBuilder::addDihedral(Molecule* molecule, vector<string>* record) {
	Compound* compound = &molecule->compounds[particle_id_maps[stoi((*record)[0])].compound_id];

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

void CompoundBuilder::countElements(Molecule* molecule) {
	int counters[4] = { 0 };
	for (int c = 0; c < molecule->n_compounds; c++) {
		Compound* C = &molecule->compounds[c];
		counters[0] += C->n_particles;
		counters[1] += C->n_pairbonds;
		counters[2] += C->n_anglebonds;
		//counters[3] += C->n_dihedrals;
	}

	printf("Molecule created with %d compounds\n", molecule->n_compounds);
	printf("%d particles added\n", counters[0]);
	printf("%d pairbonds added\n", counters[1]);
	printf("%d anglebonds added\n", counters[2]);
}





