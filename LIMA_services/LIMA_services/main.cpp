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
	Float3 pos;
	Float3 LJ_force;
	double euclidean_dist;
};

struct Row {
	Atom atoms[ATOMS_PER_ROW];
};

struct selfcenteredDatapoint {
	selfcenteredDatapoint() {}
	selfcenteredDatapoint(int id, float mass) :id(id), mass(mass){}
	int id;
	float mass;

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

void mergeSortAPI(Atom* atoms, int n_atoms) {
	Atom* sorted_atoms = mergeSort(atoms, 0, n_atoms-1);
	for (int i = 0; i < n_atoms; i++) {
		atoms[i] = sorted_atoms[i];
	}
	delete[] sorted_atoms;
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

selfcenteredDatapoint makeDatapoint(Row* rows, int atom_id, int row) {	// Row indicates row to predict, implicitly using the row before..
	Row* row_cur = &rows[row];
	Row* row_prev = &rows[row - 1];

	selfcenteredDatapoint scdp(atom_id, 0);
	for (int i = 0; i < ATOMS_PER_ROW; i++) {
		scdp.atoms_relative[i] = row_cur->atoms[i];
		scdp.atoms_relative_prev[i] = row_prev->atoms[i];
	}
	makePositionsRelative(scdp.atoms_relative, atom_id);
	makePositionsRelative(scdp.atoms_relative_prev, atom_id);

	mergeSortAPI(scdp.atoms_relative, ATOMS_PER_ROW);
	mergeSortAPI(scdp.atoms_relative_prev, ATOMS_PER_ROW);

	for (int i = 0; i < ATOMS_PER_ROW; i++) {
		//scdp.atoms_relative[i].pos.print('p');
		//scdp.atoms_relative[i].LJ_force.print('F');
	}
	

	//exit(0);
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


		row.atoms[i].pos = data[0];
		row.atoms[i].LJ_force = data[3];
		//row.atoms[i].pos.print('p');
		//row.atoms[i].LJ_force.print('F');
	}
	
	return row;
}

int readData(Row* rows) {
	string path = "D:\\Quantom\\Training\\sim_out.csv";
	fstream file;
	file.open(path);

	int row_cnt = 0;
	int lines_read = 0; 
	string line;
	int entries_per_atom = 3 * FLOAT3_PER_ATOM;	// 6 Float3







	int end_column = ATOMS_PER_ROW * 3 * FLOAT3_PER_ATOM + 10;
	string* raw_input = new string[end_column]();	// dunno about that *10
	for (int i = 0; i < end_column; i++) raw_input[i] = "";
	while (getline(file, line)) {
		int column_index = 0;
		
		if (row_cnt >= ROW_START) {

			stringstream ss(line);
			string word;
			//while (getline(ss, word, ';')) {
			//continue;
			while (getline(ss, word, ';')) {
				raw_input[column_index++] = word;
				if (column_index == end_column)
					break;
			}
			if (!((row_cnt+1) % 100))
				printf("Reading row: %d\r", row_cnt+1);
			//cout<< raw_input
			//rows[row_cnt++] = parseRow(row_raw, row_raw.size() / entries_per_atom);
			rows[lines_read++] = parseRow(raw_input);
			//printf("%d %f\n", dataline_cnt-1, rows[dataline_cnt - 1].atoms[0].LJ_force.len());
		}

		row_cnt++;

		if (row_cnt == MAX_ROW)
			break;
	}
	printf("\n");
	file.close();

	return lines_read;
}

void exportData(selfcenteredDatapoint* data, Int3 dim, string filename) {	
	std::ofstream myfile(filename);

	for (int i = 0; i < dim.x; i++) {								// Step

		if (!((i+1) % 100))
			printf("Writing row: %d\r", i+1);

		selfcenteredDatapoint scdp = data[i];
		//scdp.atoms_relative[0].LJ_force.print('f');
		scdp.atoms_relative[0].LJ_force.printToFile(&myfile);				// First print label				3xfloat
		scdp.atoms_relative_prev[0].LJ_force.printToFile(&myfile);			// Then print prev self force		3xfloat
			
		for (int ii = 1; ii < dim.y; ii++) {						
			scdp.atoms_relative[ii].pos.printToFile(&myfile);				// Print other atoms pos			69x3xfloat
			//scdp.atoms_relative[ii].LJ_force.printToFile(&myfile);
			//scdp.atoms_relative_prev[ii].pos.printToFile(&myfile);
			scdp.atoms_relative_prev[ii].LJ_force.printToFile(&myfile);		// Print other atoms prev force		69x3xfloat
		}
		if (! (i+1 == dim.x))	// Otherwise it adds an empty line at the bottom, breaking in python loading shit
			myfile << "\n";
	}
	myfile.close();
}




int main(void) {
	Row* rows = new Row[MAX_ROW + 1];
	int lines_read = readData(rows);

	selfcenteredDatapoint* data = new selfcenteredDatapoint[MAX_ROW+1];
	int cnt = 0;
	int query_atom = 0;
	printf("Processing data\n");
	for (int row = 1; row < lines_read; row += 5) {
		data[cnt++] = makeDatapoint(rows, query_atom, row);
	}

	
	exportData(data, Int3(cnt, ATOMS_PER_ROW, 1), "D:\\Quantom\\Training\\prepro_atom1.csv");



	delete[] rows;
	delete[] data;



	return 0;
}