#include "CompoundBuilder.h"





Molecule CompoundBuilder::buildMolecule(string pdb_path, string itp_path, int max_residue_id, int min_residue_id) {
	//vector<vector<string>> pdb_data = readFile(pdb_path);
	//vector<Record_ATOM> atom_data = parsePDB(pdb_path);
	compound_bridge_bundle = new CompoundBridgeBundle;
	particle_id_maps = new ParticleRef[MAX_ATOMS];

	vector<Record_ATOM> atom_data = parseGRO(pdb_path);
	vector<vector<string>> top_data = parseTOP(itp_path);


	//Compound compound;													// position of particles stored in compund.particles.pos_tsub1
	Molecule molecule;

	FFM->buildForcefield();
	//exit(0);



	//loadParticles(&compound, &atom_data, max_residue_id, true);
	loadParticles(&molecule, &atom_data, max_residue_id, min_residue_id, true);
	//printf("%d particles added\n", compound.n_particles);
	loadTopology(&molecule, &top_data);

	//molecule.compound_bridge_bundle
	molecule.compound_bridge_bundle = new CompoundBridgeBundleCompact;// (&compound_bridge_bundle);		// Convert the host version to a compact device version, belonging to the molecule
	*molecule.compound_bridge_bundle = CompoundBridgeBundleCompact(compound_bridge_bundle);

	countElements(&molecule);
	printf("Molecule built with %d compounds and %d bridges\n\n\n", molecule.n_compounds, molecule.compound_bridge_bundle->n_bridges);

	for (int i = 0; i < molecule.n_compounds; i++) {
		molecule.compounds[i].calcParticleSphere();
	}


	delete[] particle_id_maps;
	delete compound_bridge_bundle;




	return molecule;
}





void CompoundBuilder::loadParticles(Molecule* molecule, vector<CompoundBuilder::Record_ATOM>* pdb_data, int max_residue_id, int min_residue_id, bool ignore_protons) {


	int current_res_id = -1;
	int current_compound_id = -1;
	Compound* current_compound = nullptr;
	//int prev_atom_id;


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

		if (record.residue_seq_number < min_residue_id)
			continue;

		if (record.residue_seq_number > max_residue_id)
			break;

		if (ignore_protons && record.atom_name[0] == 'H')
			continue;


		



		if (record.residue_seq_number != current_res_id) {
			if (current_compound == nullptr || !current_compound->hasRoomForRes()) {
				//molecule->compound_bridge_bundle.addBridge(current_compound_id, current_compound_id + 1);
				compound_bridge_bundle->addBridge(current_compound_id, current_compound_id + 1);

				current_compound_id++;
				current_compound = &molecule->compounds[current_compound_id];
				molecule->n_compounds++;
			}
			current_res_id = record.residue_seq_number;
		}


		//particle_id_maps[record.atom_serial_number] = IDMap(record.atom_serial_number, current_compound_id, molecule->compounds[current_compound_id].n_particles);
		particle_id_maps[record.atom_serial_number] = ParticleRef(record.atom_serial_number, current_compound_id, molecule->compounds[current_compound_id].n_particles);
		//current_compound->addParticle(FFM.atomTypeToIndex(record.atom_name[0]), CompactParticle(record.coordinate));


		//ParticleParameters pp1 = FFM.getForcefield().particle_parameters[FFM.atomTypeToIndex(record.atom_name[0])];
		//ParticleParameters pp2 = FFM.getNBForcefield().particle_parameters[FFM.getAtomtypeID(record.atom_serial_number)];
		//printf("Change: %f %f %f\n\n", pp2.mass / pp1.mass, pp2.sigma / pp1.sigma, pp2.epsilon / pp1.epsilon);

		//current_compound->addParticle(FFM->atomTypeToIndex(record.atom_name[0]), record.coordinate);
		//current_compound->addParticle(FFM->getAtomtypeID(record.atom_serial_number), record.coordinate);
		current_compound->addParticle(FFM->getAtomtypeID(record.atom_serial_number), record.coordinate, FFM->atomTypeToIndex(record.atom_name[0]), record.atom_serial_number);
		molecule->n_atoms_total++;
	}
}

void CompoundBuilder::loadTopology(Molecule* molecule, vector<vector<string>>* top_data)
{
	int dihedral_cnt = 0;
	TopologyMode mode = INACTIVE;
	for (vector<string> record : *top_data) {
		if (record.size() == 0) {
			mode = INACTIVE;
			continue;
		}
		
		if (mode == INACTIVE) {
			mode = setMode(record[0]);

			if (mode == DIHEDRAL)			// Bad fix, but for now we ignore the lowest dihedral bonds, as i think they are IMPROPER DIHEDRALS
				dihedral_cnt++;

			continue;
		}
		if (mode == DIHEDRAL && dihedral_cnt > 1)
			continue;
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

void CompoundBuilder::loadMaps(ParticleRef* maps, vector<string>* record, int n) {
	for (int i = 0; i < n; i++)
		maps[i] = particle_id_maps[stoi((*record)[i])];
}

void CompoundBuilder::addGeneric(Molecule* molecule, vector<string>* record, TopologyMode mode) {
	GenericBond g_bond;
	PairBond s_bond;
	AngleBond a_bond;
	ParticleRef maps[4];

	PairBond* bondtype;
	AngleBond* angletype;
	DihedralBond* dihedraltype;


	switch (mode)
	{
	case CompoundBuilder::INACTIVE:
		break;
	case CompoundBuilder::BOND:
		addBond(molecule, maps, record);
		break;
	case CompoundBuilder::ANGLE:
		addAngle(molecule, maps, record);
		break;

	case CompoundBuilder::DIHEDRAL:
		addDihedral(molecule, maps, record);

		break;
	default:
		return;
	}
}

void CompoundBuilder::addBond(Molecule* molecule, ParticleRef* maps, vector<string>*record) {
	loadMaps(maps, record, 2);
	GenericBond g_bond = GenericBond(maps, 2);
	if (!g_bond.allParticlesExist())
		return;

	PairBond* bondtype = FFM->getBondType(maps[0].global_id, maps[1].global_id);

	distributeLJIgnores(molecule, maps, 2);

	if (!g_bond.spansTwoCompounds()) {
		Compound* compound = &molecule->compounds[maps[0].compound_id];
		compound->singlebonds[compound->n_singlebonds++] = PairBond(bondtype->b0, bondtype->kb, maps[0].local_id_compound, maps[1].local_id_compound);
		
	}
	else {
		// First, we need to make sure all bond particles are added to the bridge.
		// To create the single-bond we need to access bridge_local_indexes			
		CompoundBridge* bridge = compound_bridge_bundle->getBelongingBridge(&g_bond);
		bridge->addBondParticles(&g_bond, molecule);
		//bridge->addSinglebond(PairBond(bondtype->b0, bondtype->kb, maps[0].global_id, maps[1].global_id));
		bridge->addGenericBond(PairBond(bondtype->b0, bondtype->kb, maps[0].global_id, maps[1].global_id));
	}
}

void CompoundBuilder::addAngle(Molecule* molecule, ParticleRef* maps, vector<string>* record) {
	loadMaps(maps, record, 3);
	GenericBond g_bond = GenericBond(maps, 3);
	if (!g_bond.allParticlesExist())
		return;


	AngleBond* angletype = FFM->getAngleType(maps[0].global_id, maps[1].global_id, maps[2].global_id);

	distributeLJIgnores(molecule, maps, 3);

	if (!g_bond.spansTwoCompounds()) {
		Compound* compound = &molecule->compounds[maps[0].compound_id];
		compound->anglebonds[compound->n_anglebonds++] = AngleBond(maps[0].local_id_compound, maps[1].local_id_compound, maps[2].local_id_compound, angletype->theta_0, angletype->k_theta);
	}
	else {
		CompoundBridge* bridge = compound_bridge_bundle->getBelongingBridge(&g_bond);
		bridge->addBondParticles(&g_bond, molecule);
		bridge->addGenericBond(AngleBond(maps[0].global_id, maps[1].global_id, maps[2].global_id, angletype->theta_0, angletype->k_theta));
	}
}

void CompoundBuilder::addDihedral(Molecule* molecule, ParticleRef* maps, vector<string>* record) {
	loadMaps(maps, record, 4);

	GenericBond g_bond = GenericBond(maps, 4);
	if (!g_bond.allParticlesExist())
		return;


	DihedralBond* dihedraltype = FFM->getDihedralType(maps[0].global_id, maps[1].global_id, maps[2].global_id, maps[3].global_id);

	distributeLJIgnores(molecule, maps, 4);

	if (!g_bond.spansTwoCompounds()) {
		Compound* compound = &molecule->compounds[maps[0].compound_id];
		compound->dihedrals[compound->n_dihedrals++] = DihedralBond(maps[0].local_id_compound, maps[1].local_id_compound, maps[2].local_id_compound, maps[3].local_id_compound, dihedraltype->phi_0, dihedraltype->k_phi);
	}
	else {
		CompoundBridge* bridge = compound_bridge_bundle->getBelongingBridge(&g_bond);
		bridge->addBondParticles(&g_bond, molecule);
		bridge->addGenericBond(DihedralBond(maps[0].global_id, maps[1].global_id, maps[2].global_id, maps[3].global_id, dihedraltype->phi_0, dihedraltype->k_phi));
	}
}

void CompoundBuilder::distributeLJIgnores(Molecule* molecule, ParticleRef* particle_refs, int n) {	// Works whether particle spans multiple compounds or not.
	for (int id_self = 0; id_self < n; id_self++) {
		for (int id_other = 0; id_other < n; id_other++) {

			if (id_self == id_other)
				continue;

			Compound* compound_self = &molecule->compounds[particle_refs[id_self].compound_id];
			compound_self->lj_ignore_list[particle_refs[id_self].local_id_compound].addIgnoreTarget(particle_refs[id_other].local_id_compound, particle_refs[id_other].compound_id);
		}
	}
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
		counters[1] += C->n_singlebonds;
		counters[2] += C->n_anglebonds;
		counters[3] += C->n_dihedrals;
	}

	printf("Molecule created with %d compounds\n", molecule->n_compounds);
	printf("%d particles added\n", counters[0]);
	printf("%d singlebonds added\n", counters[1]);
	printf("%d anglebonds added\n", counters[2]);
	printf("%d dihedrals added\n", counters[3]);
}















Molecule::Molecule() {
	compounds = new Compound[MAX_COMPOUNDS];
	//compound_bridge_bundle = new CompoundBridgeBundleCompact;
}


