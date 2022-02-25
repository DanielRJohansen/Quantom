#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;


const int ATOMS_PER_ROW = 70;	// How many we are interested in loading. This loading is BEFORE we sort by distance
const int FLOAT3_PER_ATOM = 6;
const int ROW_START = 100;
const int MAX_ROW = 1e+5;
//const int MAX_ROW = 1e+3;


struct Int3 {
	Int3(){}
	Int3(int x, int y, int z) : x(x), y(y), z(z) {}

	int x, y, z;
};

struct Float3 {
	Float3() {}
	Float3(double a) : x(a), y(a), z(a) {}
	Float3(double x, double y, double z) : x(x), y(y), z(z) {}
	Float3(double* a) { x = a[0]; y = a[1]; z = a[2]; }

	inline Float3 operator * (const double a) const { return Float3(x * a, y * a, z * a); }
	inline Float3 operator * (const Float3 a) const { return Float3(x * a.x, y * a.y, z * a.z); }
	inline Float3 operator + (const Float3 a) const { return Float3(x + a.x, y + a.y, z + a.z); }
	inline Float3 operator - (const Float3 a) const { return Float3(x - a.x, y - a.y, z - a.z); }
	inline bool operator == (const Float3 a) const { return (a.x == x && a.y == y && a.z == z); }
	inline void operator += (const Float3 a) { x += a.x; y += a.y; z += a.z; }
	inline void operator *= (const double a) { x *= a; y *= a; z *= a; }

	inline bool operator < (const Float3 a) { return x < a.x&& y < a.y&& z < a.z; }
	inline bool operator > (const Float3 a) { return x > a.x && y > a.y && z > a.z; }

	void print(char c = '_') { printf("%c %f %f %f\n", c, x, y, z); }


	Float3 norm() {
		double l = len();
		if (l)
			return *this * (1.f / l);
		return Float3(0, 0, 0);
	}
	Float3 square() { return Float3(x * x, y * y, z * z); }
	inline double len() { return (double)sqrtf(x * x + y * y + z * z); }
	double at(int index) {
		switch (index) {
		case 0:
			return x;
		case 1:
			return y;
		case 2:
			return z;
		default:

			return -404;
		}
	}
	double* placeAt(int index) {
		switch (index) {
		case 0:
			return &x;
		case 1:
			return &y;
		case 2:
			return &z;
		}
	}

	void printToFile(std::ofstream* open_file) {
		//print('s');
		for (int j = 0; j < 3; j++) {
			*open_file << this->at(j) << ';';
		}
	}

	double x, y, z;
};

struct Atom {
	int id;
	Float3 pos;
	Float3 LJ_force;
	double euclidean_dist;
};

struct Row {
	Atom atoms[ATOMS_PER_ROW];
};

struct selfcenteredDatapoint {
	selfcenteredDatapoint() {}
	selfcenteredDatapoint(int id, float mass) :queryatom_id(id), mass(mass){}
	int queryatom_id;
	float mass;

	bool ignore = false;

	Atom atoms_relative[128];	// Sorted with closest neighbor first
	Atom atoms_relative_prev[128];	// Sorted with closest neighbor first
};


Atom* merge(const Atom* left, int n_left, const Atom* right, int n_right) {
	Atom* merged = new Atom[n_left + n_right];
	int l = 0;
	int r = 0;
	int index = 0;

	while (l < n_left && r < n_right) {
		if (left[l].euclidean_dist < right[r].euclidean_dist) {
			merged[index++] = left[l++];
		}
		else {
			merged[index++] = right[r++];
		}
	}
	while (l < n_left) { merged[index++] = left[l++]; }
	while (r < n_right) { merged[index++] = right[r++]; }


	return merged;
}

Atom* mergeSort(const Atom* atoms, const int l, const int r) {	// l and r are indexes of the two extremes
	int m = (l + r) / 2;

	Atom* left;
	Atom* right;

	//printf("step\n");
	if (r - l > 1) {
		//printf("l %d r %d\n", l, r);
		left = mergeSort(atoms, l, m-1);
		right = mergeSort(atoms, m, r);
		Atom* merged = merge(left, m - l, right, r - m + 1);

		if (l == m - 1)		// Take care of special case, where only left side can be a single object instead of array!
			delete left;
		else
			delete[] left;
		delete[] right;

		return merged;
	}
	else if (r - l == 1) {
		return merge(&atoms[l], 1, &atoms[r], 1);
	}
	else {
		return new Atom(atoms[l]);
	}
}

int* mergeSortAPI(Atom* atoms, int n_atoms) {					// Returns a mapping where mapping[0] is the closest id, mapping [1] seconds closest, so on
	Atom* sorted_atoms = mergeSort(atoms, 0, n_atoms-1);
	for (int i = 0; i < n_atoms; i++) {
		atoms[i] = sorted_atoms[i];
	}
	
	int* mapping = new int[n_atoms];
	for (int i = 0; i < n_atoms; i++) mapping[i] = sorted_atoms[i].id;

	delete[] sorted_atoms;
	return mapping;
}

void sortAtoms(Atom* atoms, int n_atoms, int* mapping) {
	Atom* temp = new Atom[n_atoms];
	for (int i = 0; i < n_atoms; i++) {
		temp[i] = atoms[i];
	}
	for (int i = 0; i < n_atoms; i++) {
		temp[i] = atoms[mapping[i]];
	}
	delete[] temp;
}


void makePositionsRelative(Atom* atoms, int atom_id) {
	for (int i = 0; i < ATOMS_PER_ROW; i++) {
		if (i == atom_id)
			continue;
		atoms[i].pos = atoms[i].pos - atoms[atom_id].pos;
		atoms[i].euclidean_dist = atoms[i].pos.len();
	}
	atoms[atom_id].euclidean_dist = 0;	// The heliocentric atom always comes first
	atoms[atom_id].pos = Float3(0.f);
}

selfcenteredDatapoint makeDatapoint(Row* rows, int query_atom_id, int row) {	// Row indicates row to predict, implicitly using the row before..
	Row* row_cur = &rows[row];
	Row* row_prev = &rows[row - 1];

	selfcenteredDatapoint scdp(query_atom_id, 0);
	for (int i = 0; i < ATOMS_PER_ROW; i++) {			// Load positions. They are NOT yet relative to query atom
		scdp.atoms_relative[i] = row_cur->atoms[i];
		scdp.atoms_relative_prev[i] = row_prev->atoms[i];
	}
	makePositionsRelative(scdp.atoms_relative, query_atom_id);
	makePositionsRelative(scdp.atoms_relative_prev, query_atom_id);

	int* mapping = mergeSortAPI(scdp.atoms_relative, ATOMS_PER_ROW);
	sortAtoms(scdp.atoms_relative, ATOMS_PER_ROW, mapping);
	sortAtoms(scdp.atoms_relative_prev, ATOMS_PER_ROW, mapping);
	delete[] mapping;

	
	return scdp;
}



Row parseRow(string* raw_data) {
	int index = 0;
	Row row;
	Float3 data[6];
	//return Row();
	for (int i = 0; i < ATOMS_PER_ROW; i++) {
		for (int ii = 0; ii < FLOAT3_PER_ATOM; ii++) {
			for (int iii = 0; iii < 3; iii++) {
				//cout << raw_data[iii + ii * 3 + i * FLOAT3_PER_ATOM * 3] << endl;
				*data[ii].placeAt(iii) = stod(raw_data[iii + ii * 3 + i * FLOAT3_PER_ATOM * 3]);
			}
		}

		row.atoms[i].id = i;
		row.atoms[i].pos = data[0];
		row.atoms[i].LJ_force = data[3];
		//row.atoms[i].pos.print('p');
		//row.atoms[i].LJ_force.print('F');
	}
	
	return row;
}


void readDataBIN(Row* rows, string path, int N_STEPS) {
	char* file_path;
	file_path = &path[0];
	cout << "Reading from file " << file_path << endl;

	FILE* file;
	fopen_s(&file, file_path, "rb");

	int particles = 128;
	int f3_per_particle = 6;
	int n_f3s = N_STEPS * particles * f3_per_particle;
	Float3* buffer = new Float3[n_f3s];
	fread(buffer, sizeof(Float3), n_f3s, file);
	fclose(file);


	for (int i = 0; i < N_STEPS; i++) {
		for (int j = 0; j < particles; j++) {
			if (j >= 70)
				continue;
			int buffer_index = j * f3_per_particle + i * particles * f3_per_particle;
			rows[i].atoms[j].pos = buffer[buffer_index + 0];
			rows[i].atoms[j].LJ_force = buffer[buffer_index + 3];
			rows[i].atoms[j].id = j;
		}
	}
}

void exportDataCSV(selfcenteredDatapoint* data, Int3 dim, string filename) {
	std::ofstream myfile(filename);
	printf("\nExporting to file: \t");
	cout << filename << endl;



	for (int i = 0; i < dim.x; i++) {								// Step
		if ((i % 100) == 99)
			printf("\rWriting row: %d", i + 1);

		selfcenteredDatapoint scdp = data[i];
		if (scdp.ignore)
			continue;

		//scdp.atoms_relative[0].LJ_force.print('f');
		scdp.atoms_relative[0].LJ_force.printToFile(&myfile);				// First print label				3xfloat
		scdp.atoms_relative_prev[0].LJ_force.printToFile(&myfile);			// Then print prev self force		3xfloat

		for (int ii = 1; ii < dim.y; ii++) {
			scdp.atoms_relative[ii].pos.printToFile(&myfile);				// Print other atoms pos			69x3xfloat
			//scdp.atoms_relative[ii].LJ_force.printToFile(&myfile);
			//scdp.atoms_relative_prev[ii].pos.printToFile(&myfile);
			scdp.atoms_relative_prev[ii].LJ_force.printToFile(&myfile);		// Print other atoms prev force		69x3xfloat
		}
		if (!(i + 1 == dim.x))	// Otherwise it adds an empty line at the bottom, breaking in python loading shit
			myfile << "\n";
	}
	myfile.close();
}

void exportData(selfcenteredDatapoint* data, Int3 dim, string path) {				// 	exportData(data, Int3(n_datapoints, ATOMS_PER_ROW, 1), path_out);
	char* file_path;
	file_path = &path[0];
	cout << "Writing to file " << file_path << endl;

	int lines_to_print = 0;	
	Float3* buffer = new Float3[dim.x * dim.y * 10];	// Times 10 just to be safe
	uint32_t buffer_ptr = 0;

	for (int i = 0; i < dim.x; i++) {								// Step
		selfcenteredDatapoint scdp = data[i];
		if (scdp.ignore)
			continue;


		//scdp.atoms_relative[0].LJ_force.print('f');
		buffer[buffer_ptr++] = scdp.atoms_relative[0].LJ_force;					// First print label				3xfloat
		buffer[buffer_ptr++] = scdp.atoms_relative_prev[0].LJ_force;			// Then print prev self force		3xfloat
			
		for (int ii = 1; ii < dim.y; ii++) {						
			buffer[buffer_ptr++] = scdp.atoms_relative[ii].pos;				// Print other atoms pos			69x3xfloat
			buffer[buffer_ptr++] = scdp.atoms_relative_prev[ii].LJ_force;		// Print other atoms prev force		69x3xfloat
		}
	}
	
	printf("Ptr val: %u\n", buffer_ptr);

	FILE* file;
	fopen_s(&file, file_path, "wb");
	fwrite(buffer, sizeof(Float3), buffer_ptr, file);
	fclose(file);
}

void makeForceChangePlot(selfcenteredDatapoint* data, int n_datapoints) {
	int bins[32] = { 0 };
	int zeroes = 0;
	int bin;
	for (int i = 0; i < n_datapoints; i++) {
		double force_change = (data[i].atoms_relative[0].LJ_force - data[i].atoms_relative_prev[0].LJ_force).len();
		
		if (force_change < 1.f) {
			bin = 0;
			if (force_change == 0.f) {
				zeroes++;
				data[i].atoms_relative[1].pos.x = 40404;
				data[i].atoms_relative[0].LJ_force = Float3(40404.f);
				printf("What");
			}
				
		}
		else 
			bin = (int)log2(force_change);
			
		//printf("BIN: %d. Force %f\n", bin, force_change);
		bins[bin]++;
		if (bin < 0 || bin > 31) {
			printf("Bin: %d\n");
			exit(0);
		}
			
	}
	printf("\n%d zeroes in dataset\n", zeroes);

	printf("x=categorical([");
	for (int i = 0; i < 32; i++) printf("%u ", (int) pow(2, i));
	printf("]);\ny=[");
	for (int i = 0; i < 32; i++) printf("%d ", bins[i]);
	printf("];\n\n\n");
}

int discardVolatileDatapoints(selfcenteredDatapoint* data, int n_datapoints, double cuttoff_value) {	// Removes datapoints where the change in force is too large for the NN to handle..
	int n_ignored = 0;
	for (int i = 0; i < n_datapoints; i++) {
		double force_change = (data[i].atoms_relative[0].LJ_force - data[i].atoms_relative_prev[0].LJ_force).len();
		if (force_change > cuttoff_value) {
			data[i].ignore = true;
			n_ignored++;
		}
	}
	printf("\n%d datapoints cuttoff due to force change magnitude\n", n_ignored);
	return n_ignored;
}


#include <algorithm>

void shuffle(selfcenteredDatapoint* data, int n_datapoints) {
	int* order = new int[n_datapoints];
	for (int i = 0; i < n_datapoints; i++) order[i] = i;
	random_shuffle(order, order + n_datapoints);
	

	selfcenteredDatapoint* temp = new selfcenteredDatapoint[n_datapoints];
	for (int i = 0; i < n_datapoints; i++) temp[i] = data[i];
	for (int i = 0; i < n_datapoints; i++) data[i] = temp[order[i]];
	delete[] order;
}




int main(int argc, char** argv) {


	//cout << argv[0] << endl << argv[1] << endl << argv[2] << endl << "DONE" << endl;
	printf("argc %d\n", argc);
	for (int i = 0; i < argc; i++) {
		cout << ":" << argv[i] << ":" << endl;
	}

	string workdir = "D:\\Quantom\LIMANET\sim_out";
	bool shuffle_time_dim = false;
	int N_STEPS = 100000;	// Determines file to read

	printf("argc %d\n", argc);
	if (argc > 0) {
		workdir = argv[1];
		N_STEPS = stoi(argv[2]);
		shuffle_time_dim = stoi(argv[3]);
	}
		



	printf("Allocating %.01f GB of RAM\n", ((double)sizeof(Row) * MAX_ROW + (double) sizeof(selfcenteredDatapoint) * MAX_ROW)/ 1000000000.f);
	Row* rows = new Row[MAX_ROW + 1];
	selfcenteredDatapoint* data = new selfcenteredDatapoint[MAX_ROW + 1];
	int n_datapoints = 0;
	int query_atom = 0;
	



	printf("Processing data\n");



	string path_in = workdir + "\\sim_out.bin";
	readDataBIN(rows, path_in, N_STEPS);
	for (int row = ROW_START; row < N_STEPS; row += 2) {
		data[n_datapoints++] = makeDatapoint(rows, query_atom, row);
	}




	makeForceChangePlot(data, n_datapoints);
	int n_ignored_datapoints = discardVolatileDatapoints(data, n_datapoints, 5000.f);


	//string path_out = "D:\\Quantom\\LIMANET\\sim_out\\atom" + to_string(query_atom) + "_lines" + to_string(n_datapoints - n_ignored_datapoints);
	string path_out = workdir + "\\traindata";
	if (shuffle_time_dim) {
		shuffle(data, n_datapoints);
		path_out = path_out + "_shuffled";
	}
	exportData(data, Int3(n_datapoints, ATOMS_PER_ROW, 1), path_out + ".bin");



	delete[] rows;
	delete[] data;

	return 0;
}