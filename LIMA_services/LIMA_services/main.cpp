#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;


const int FLOAT3_PER_ATOM = 3;
const int ROW_START = 500;
const int MAX_NEIGHBORS_OUT = 32;		// AFter sorting, only nearest x neighbors is printed to the out-file
//const int MAX_ROW = 1e+3;
float BOX_LEN = 7.f;
float BOX_LEN_HALF = BOX_LEN/2.f;




float max_dif = 0;

struct Float3 {
	Float3() {}
	Float3(float a) : x(a), y(a), z(a) {}
	Float3(float x, float y, float z) : x(x), y(y), z(z) {}
	Float3(float* a) { x = a[0]; y = a[1]; z = a[2]; }

	inline Float3 operator * (const float a) const { return Float3(x * a, y * a, z * a); }
	inline Float3 operator * (const Float3 a) const { return Float3(x * a.x, y * a.y, z * a.z); }
	inline Float3 operator + (const Float3 a) const { return Float3(x + a.x, y + a.y, z + a.z); }
	inline Float3 operator - (const Float3 a) const { return Float3(x - a.x, y - a.y, z - a.z); }
	inline bool operator == (const Float3 a) const { return (a.x == x && a.y == y && a.z == z); }
	inline void operator += (const Float3 a) { x += a.x; y += a.y; z += a.z; }
	inline void operator *= (const float a) { x *= a; y *= a; z *= a; }

	inline bool operator < (const Float3 a) { return x < a.x&& y < a.y&& z < a.z; }
	inline bool operator > (const Float3 a) { return x > a.x && y > a.y && z > a.z; }

	void print(char c = '_') { printf("%c %f %f %f\n", c, x, y, z); }


	Float3 norm() {
		float l = len();
		if (l)
			return *this * (1.f / l);
		return Float3(0, 0, 0);
	}
	Float3 square() { return Float3(x * x, y * y, z * z); }
	inline float len() { return (float)sqrtf(x * x + y * y + z * z); }
	float at(int index) {
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
	float* placeAt(int index) {
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

	float _max() {
		return max(x, max(y, z));
	}
	float _min() {
		return min(x, max(y, z));
	}

	float x, y, z;
	//float x=0, y=0, z=0;
};

struct Atom {
	Atom(){}
	Atom(int id, Float3 pos, Float3 LJ_force) : id(id), pos(pos), LJ_force(LJ_force){}
	int id;
	Float3 pos;
	Float3 LJ_force;
	float euclidean_dist;
};

struct Row {
	Row(){}
	Row(int atoms_per_row) {
		atoms = new Atom[atoms_per_row];
	}

	static int byteSize(int atoms_per_row) {
		return sizeof(Atom) * atoms_per_row;
	}

	~Row() {
		//delete[] atoms;
	}
	Atom* atoms;
};

struct selfcenteredDatapoint {
	selfcenteredDatapoint() {}
	selfcenteredDatapoint(int id,int particles_per_step) :queryatom_id(id){
		atoms_relative = new Atom[particles_per_step];
		atoms_relative_tsub1 = new Atom[particles_per_step];
		atoms_relative_tsub2 = new Atom[particles_per_step];
		initialized = true;
	}
	int queryatom_id;

	bool ignore = false;
	bool initialized = false;


	Atom* atoms_relative;	// Sorted with closest neighbor first
	Atom* atoms_relative_tsub1;	// Sorted with closest neighbor first
	Atom* atoms_relative_tsub2;	// Sorted with closest neighbor first

	static int byteSize(int particles_per_step) {
		return sizeof(Atom) * particles_per_step * 3 + sizeof(selfcenteredDatapoint);
	}

	void freeMem() {
		delete[] atoms_relative;
		delete[] atoms_relative_tsub1;
		delete[] atoms_relative_tsub2;
	}


};

void applyHyperpos(Float3* static_particle, Float3* movable_particle) {
	for (int i = 0; i < 3; i++) {
		*movable_particle->placeAt(i) += BOX_LEN * ((static_particle->at(i) - movable_particle->at(i)) > BOX_LEN_HALF);
		*movable_particle->placeAt(i) -= BOX_LEN * ((static_particle->at(i) - movable_particle->at(i)) < -BOX_LEN_HALF);	// use at not X!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	}
}

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


void makePositionsRelative(Atom* atoms, int atom_id, int particles_per_step) {
	for (int i = 0; i < particles_per_step; i++) {
		if (i == atom_id)
			continue;

		applyHyperpos(&atoms[atom_id].pos, &atoms[i].pos);
		atoms[i].pos = atoms[i].pos - atoms[atom_id].pos;
		atoms[i].euclidean_dist = atoms[i].pos.len();
	}
	atoms[atom_id].euclidean_dist = 0;	// The heliocentric atom always comes first
	atoms[atom_id].pos = Float3(0.f);	// If we erase this pos, we can't get delta positions
}

selfcenteredDatapoint makeDatapoint(Row* rows, int query_atom_id, int row, int particles_per_step) {	// Row indicates row to predict, implicitly using the row before..
	Row* row_cur = &rows[row];
	Row* row_tsub1 = &rows[row - 1];
	Row* row_tsub2 = &rows[row - 2];

	selfcenteredDatapoint scdp(query_atom_id, particles_per_step);
	for (int i = 0; i < particles_per_step; i++) {			// Load positions. They are NOT yet relative to query atom
		scdp.atoms_relative[i] = row_cur->atoms[i];
		scdp.atoms_relative_tsub1[i] = row_tsub1->atoms[i];
		scdp.atoms_relative_tsub2[i] = row_tsub2->atoms[i];
	}
	makePositionsRelative(scdp.atoms_relative, query_atom_id, particles_per_step);
	makePositionsRelative(scdp.atoms_relative_tsub1, query_atom_id, particles_per_step);
	makePositionsRelative(scdp.atoms_relative_tsub2, query_atom_id, particles_per_step);

	int* mapping = mergeSortAPI(scdp.atoms_relative, particles_per_step);		// First generation a sorting mapping
	sortAtoms(scdp.atoms_relative, particles_per_step, mapping);						// Use mapping to sort both data-arrays, so things match!
	sortAtoms(scdp.atoms_relative_tsub1, particles_per_step, mapping);
	sortAtoms(scdp.atoms_relative_tsub2, particles_per_step, mapping);
	delete[] mapping;

	
	return scdp;
}


void readDataBIN(Row* rows, string path, int N_STEPS, int particles_per_step) {
	char* file_path;
	file_path = &path[0];
	cout << "Reading from file " << file_path << endl;

	FILE* file;
	fopen_s(&file, file_path, "rb");

	
	uint64_t n_f3s = N_STEPS * particles_per_step * FLOAT3_PER_ATOM;
	printf("Allocating data-buffer of size %.2Lf GB for reading\n", (long double)n_f3s * sizeof(Float3) * 1e-9);
	Float3* buffer = new Float3[n_f3s];
	fread(buffer, sizeof(Float3), n_f3s, file);
	fclose(file);


	for (int step = 0; step < N_STEPS; step++) {
		if (!(step % 1000))
			printf("\rReading row %d", step);

		for (int p = 0; p < particles_per_step; p++) {			
			uint64_t buffer_index = p * FLOAT3_PER_ATOM + step * particles_per_step * FLOAT3_PER_ATOM;

			rows[step].atoms[p] = Atom(p, buffer[buffer_index + 0], buffer[buffer_index + 1]);	// Need of way of determining whether a particle exists, or is in the buffer!!!!!!!!!!!!!!!!!!!!!
		}
	}

	delete[] buffer;
	printf("\n");
}



void appendToBuffer(selfcenteredDatapoint* scdp, int particles_per_datapoint, int max_neighbors, Float3* buffer) {
	//int lines_to_print = 0;
	//Float3* buffer = new Float3[n_datapoints * 2 * (1 + particles_per_datapoint)];	// Times 10 just to be safe
	
	uint32_t buffer_ptr_rel = 0;

	if (scdp->ignore)	
		return;

	float force_dif = scdp->atoms_relative[0].LJ_force.x - scdp->atoms_relative_tsub1[0].LJ_force.x;
	max_dif = max(abs(force_dif), max_dif);
	


	//scdp->atoms_relative[0].LJ_force.print('f');
	buffer[buffer_ptr_rel++] = scdp->atoms_relative[0].LJ_force;													// First print label				
	buffer[buffer_ptr_rel++] = scdp->atoms_relative[0].pos - scdp->atoms_relative_tsub1[0].pos;					// Then vel
	buffer[buffer_ptr_rel++] = scdp->atoms_relative_tsub1[0].LJ_force;											// Then print prev self force		
	buffer[buffer_ptr_rel++] = scdp->atoms_relative_tsub1[0].LJ_force - scdp->atoms_relative_tsub2[0].LJ_force;	// Force change



	for (int ii = 1; ii < particles_per_datapoint; ii++) {
		buffer[buffer_ptr_rel++] = scdp->atoms_relative[ii].pos;				// Print other atoms pos			69x3xfloat
		buffer[buffer_ptr_rel++] = (scdp->atoms_relative[ii].pos - scdp->atoms_relative_tsub2[ii].pos) * 0.5;
		buffer[buffer_ptr_rel++] = scdp->atoms_relative_tsub1[ii].LJ_force;		// Print other atoms prev force		69x3xfloat
		buffer[buffer_ptr_rel++] = scdp->atoms_relative_tsub1[ii].LJ_force - scdp->atoms_relative_tsub2[ii].LJ_force;

		if (ii == max_neighbors)
			break;
	}
	

	//printf("\nPtr val: %u\n", buffer_ptr);

}




void exportData( string path, Float3* buffer, uint64_t buffer_len) {				// 	exportData(data, Int3(n_datapoints, ATOMS_PER_ROW, 1), path_out);
	char* file_path;
	file_path = &path[0];
	cout << "Writing to file " << file_path << endl;

	/*
	int lines_to_print = 0;	
	Float3* buffer = new Float3[n_datapoints * 2 * (1 + particles_per_datapoint)];	// Times 10 just to be safe
	uint32_t buffer_ptr = 0;

	for (int step = 0; step < n_datapoints; step++) {								// Step
		selfcenteredDatapoint scdp = data[step];
		if (scdp.ignore)
			continue;


		//scdp.atoms_relative[0].LJ_force.print('f');
		buffer[buffer_ptr++] = scdp.atoms_relative[0].LJ_force;													// First print label				
		buffer[buffer_ptr++] = scdp.atoms_relative[0].pos - scdp.atoms_relative_tsub1[0].pos;					// Then vel
		buffer[buffer_ptr++] = scdp.atoms_relative_tsub1[0].LJ_force;											// Then print prev self force		
		buffer[buffer_ptr++] = scdp.atoms_relative_tsub1[0].LJ_force - scdp.atoms_relative_tsub2[0].LJ_force;	// Force change
		
	

		for (int ii = 1; ii < particles_per_datapoint; ii++) {
			buffer[buffer_ptr++] = scdp.atoms_relative[ii].pos;				// Print other atoms pos			69x3xfloat
			buffer[buffer_ptr++] = (scdp.atoms_relative[ii].pos - scdp.atoms_relative_tsub2[ii].pos) * 0.5;
			buffer[buffer_ptr++] = scdp.atoms_relative_tsub1[ii].LJ_force;		// Print other atoms prev force		69x3xfloat
			buffer[buffer_ptr++] = scdp.atoms_relative_tsub1[ii].LJ_force - scdp.atoms_relative_tsub2[ii].LJ_force;

			if (ii == max_neighbors)
				break;
		}
	}
	
	printf("Ptr val: %u\n", buffer_ptr);


	float maxVal = 0; 
	float minVal = 0; 
	for (int i = 0; i < buffer_ptr; i++) {
		maxVal = max(maxVal, buffer[i]._max());
		minVal = min(minVal, buffer[i]._min());
	}

	printf("Max: %f Min: %f\n", maxVal, minVal);
	*/


	FILE* file;
	fopen_s(&file, file_path, "wb");
	//fwrite(buffer, sizeof(Float3), buffer_ptr, file);
	fwrite(buffer, sizeof(Float3), buffer_len, file);
	fclose(file);
}

void makeForceChangePlot(selfcenteredDatapoint* data, int n_datapoints) {
	const int n_bins = 22;

	int bins[n_bins] = { 0 };
	int zeroes = 0;
	int bin;
	for (int i = 0; i < n_datapoints; i++) {
		float force_change = (data[i].atoms_relative[0].LJ_force - data[i].atoms_relative_tsub1[0].LJ_force).len();				// Only the forces changes that the NN must predict!
		
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
	for (int i = 0; i < n_bins; i++) printf("%u ", (int) pow(2, i));
	printf("]);\ny=[");
	for (int i = 0; i < n_bins; i++) printf("%d ", bins[i]);
	printf("];\n\n\n");
}

int discardVolatileDatapoints(selfcenteredDatapoint* data, int n_datapoints, float cuttoff_value) {	// Removes datapoints where the change in force is too large for the NN to handle..
	int n_ignored = 0;
	for (int i = 0; i < n_datapoints; i++) {
		float force_change = (data[i].atoms_relative[0].LJ_force - data[i].atoms_relative_tsub1[0].LJ_force).len();		// It's okay we only do with the one that the NN has to predict!
		if (force_change > cuttoff_value) {
			data[i].ignore = true;
			n_ignored++;
		}
	}
	printf("\n%d datapoints cuttoff due to force change magnitude\n", n_ignored);
	return n_ignored;
}

int discardVolatileDatapoint(selfcenteredDatapoint* scdp, float cuttoff_value) {	// Removes datapoints where the change in force is too large for the NN to handle..
	float force_change = (scdp->atoms_relative[0].LJ_force - scdp->atoms_relative_tsub1[0].LJ_force).len();		// It's okay we only do with the one that the NN has to predict!
	if (force_change > cuttoff_value) {
		scdp->ignore = true;
		return 1;
	}
	
	//printf("\n%d datapoints cuttoff due to force change magnitude\n", n_ignored);
	return 0;
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
	delete[] temp;
}



// This code should only be run alone for testing! Always let Quantom call this program!!
int main(int argc, char** argv) {


	//cout << argv[0] << endl << argv[1] << endl << argv[2] << endl << "DONE" << endl;
	printf("argc %d\n", argc);
	for (int i = 0; i < argc; i++) {
		cout << ":" << argv[i] << ":" << endl;
	}

	//string workdir = "D:\\Quantom\LIMANET\sim_out";
	int N_STEPS = 350000;	// Determines file to read
	string workdir = "C:\\PROJECTS\\Quantom\\Simulation\\Steps_" + to_string(N_STEPS) + "\\";
	
	bool shuffle_time_dim = false;
	int n_compounds = 13;
	int particles_per_compound = 128;



	printf("argc %d\n", argc);
	if (argc > 5) {
		workdir = argv[1];
		N_STEPS = stoi(argv[2]);
		shuffle_time_dim = stoi(argv[3]);
		n_compounds = stoi(argv[4]);
		int particles_per_compound = stoi(argv[5]);
	}
	else if (argc != 1) {
		printf("Not enough input arguments for processing\n");
		exit(1);
	}
		
	int particles_per_step = particles_per_compound * n_compounds;





	printf("Allocating %.02Lf GB of RAM for rows\n", (long double) Row::byteSize(particles_per_step) * N_STEPS * 1e-9 );	// Gives approx size 3 = 1 in row, 2 in scdp
	Row* rows = new Row[N_STEPS];
	for (int i = 0; i < N_STEPS; i++)
		rows[i] = Row(particles_per_step);

	//selfcenteredDatapoint* data = new selfcenteredDatapoint[N_STEPS];
	int n_datapoints = 0;
	int query_atom = 0;
	



	printf("Processing data\n");



	string path_in = workdir + "sim_traindata.bin";
	readDataBIN(rows, path_in, N_STEPS, particles_per_step);




	int spacing = 2;
	int n_queries = 2;
	//uint64_t datapoints_per_query = ceil((double)N_STEPS - (double)ROW_START) / (double)spacing;
	//uint64_t total_datapoints = datapoints_per_query * n_queries;
	int f3_per_datapoint = 4 + 4 * MAX_NEIGHBORS_OUT;
	uint64_t steps_total_out = ceil((double)N_STEPS - (double)ROW_START) / (double)spacing;
	uint64_t f3_per_step = f3_per_datapoint * n_queries;
	uint64_t f3_total = steps_total_out * f3_per_step;


	//uint64_t f3_per_query = f3_per_datapoint * datapoints_per_query;


	//printf("Allocating %.02Lf GB of RAM for out-buffer\n", (long double)sizeof(Float3) * total_datapoints * f3_per_datapoint * 1e-9);	// Gives approx size 3 = 1 in row, 2 in scdp
	printf("Allocating %.02Lf GB of RAM for out-buffer\n", (long double)sizeof(Float3) * steps_total_out * f3_per_step * 1e-9);	// Gives approx size 3 = 1 in row, 2 in scdp
	Float3* buffer = new Float3[steps_total_out * f3_per_step];	// Times 10 just to be safe

	int step_index = 0;
	for (int row = ROW_START; row < N_STEPS; row += spacing) {
		if (!((row - ROW_START) % 100))
			printf("\rProcessing datapoint at row %07d", row);

		for (int query_id = 0; query_id < n_queries; query_id++) {
			
			selfcenteredDatapoint scdp = makeDatapoint(rows, query_id, row, particles_per_step);

			//uint64_t buffer_index = query_id * f3_per_query + step_index * f3_per_datapoint;

			uint64_t buffer_index = query_id * f3_per_datapoint + step_index * f3_per_step;
			discardVolatileDatapoint(&scdp, 1e+4);
			//printf("Buffer_in")
			appendToBuffer(&scdp, particles_per_step, MAX_NEIGHBORS_OUT, &buffer[buffer_index]);
			scdp.freeMem();
		}
		step_index++;
		//printf("\n Pos of n1 %f\n", data[0].atoms_relative[3].pos.x);

		//appendToBuffer(data, datapoints_per_query, particles_per_step, MAX_NEIGHBORS_OUT, &buffer[query_id * f3_per_query]);
	}
	printf("\n\n");

	printf("Max force dif %f\n", max_dif);

	string path_out = workdir + "\\traindata_queries" + to_string(n_queries);

	exportData(path_out + ".bin", buffer, f3_total);



	/*


	for (int row = ROW_START; row < N_STEPS; row += spacing) {
		if (!((row-ROW_START) % 100))
			printf("\rProcessing datapoint at row %07d", row);
		data[n_datapoints++] = makeDatapoint(rows, query_atom, row, particles_per_step);
	}
	printf("\n%d datapoints \n\n", n_datapoints);



	//makeForceChangePlot(data, n_datapoints);
	//int n_ignored_datapoints = discardVolatileDatapoints(data, n_datapoints, 1000.f);


	//string path_out = "D:\\Quantom\\LIMANET\\sim_out\\atom" + to_string(query_atom) + "_lines" + to_string(n_datapoints - n_ignored_datapoints);
	string path_out = workdir + "\\traindata1";
	if (shuffle_time_dim) {
		shuffle(data, n_datapoints);
		path_out = path_out + "_shuffled";
	}

	uint64_t buffer_ptr_total = 0;
	Float3* buffer = new Float3[n_datapoints * 2 * (1 + particles_per_step)];	// Times 10 just to be safe
	appendToBuffer(data, n_datapoints, particles_per_step, MAX_NEIGHBORS_OUT, buffer, &buffer_ptr_total);
	exportData(data, n_datapoints, particles_per_step, path_out + ".bin", MAX_NEIGHBORS_OUT, buffer, buffer_ptr_total);
	*/


	delete[] rows;			// We are NOT properly delete[]'ing atoms in row
	delete[] buffer;
	return 0;
}