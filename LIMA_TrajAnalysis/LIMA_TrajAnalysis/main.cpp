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
float BOX_LEN_HALF = BOX_LEN / 2.f;

const float PI = 3.1415;


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
	float dot(Float3 a) const { return (x * a.x + y * a.y + z * a.z); }
	Float3 cross(Float3 a) const { return Float3(y * a.z - z * a.y, z * a.x - x * a.z, x * a.y - y * a.x); }
	static float getAngle(Float3 v1, Float3 v2) {
		float val = (v1.dot(v2)) / (v1.len() * v2.len());	// If i make this float, we get values over 1, even with the statements below! :(
		//if (val > 1.f || val < -1.f) { printf("Val1 %f !!\n", val);}
		val = val > 1.f ? 1.f : val;
		val = val < -1.f ? -1.f : val;
		if (val > 1.f || val < -1.f) {
			printf("Val2 %f !!\n", val);
		}
		return acos(val);
	}
	static Float3 rodriguesRotatation(Float3 v, Float3 k, double theta) {
		return v * cos(theta) + k.cross(v) * sin(theta) + k * (k.dot(v)) * (1 - cos(theta));
	}

	float x, y, z;
	//float x=0, y=0, z=0;
};

struct Particle {
	Particle() {}
	Particle(int id, Float3 pos, Float3 LJ_force) : id(id), pos(pos), LJ_force(LJ_force) {}
	int id;
	Float3 pos;
	Float3 LJ_force;
	float euclidean_dist;
	bool inactive = false;
};

struct Row {
	Row() {}
	Row(int particles_per_row) {		
		n_particles = particles_per_row;
		particles = new Particle[n_particles];
	}

	static int byteSize(int particles_per_row) {
		return sizeof(Particle) * particles_per_row;
	}

	~Row() {
		//delete[] particles;
	}
	Particle* particles;
	int n_particles;
};

Float3 calcCenterOfMass(Row* row) {
	Float3 sum = Float3(0.f);
	int n_actives = 0;

	for (int i = 0; i < row->n_particles; i++) {
		if (!row->particles[i].inactive) {
			sum += row->particles[i].pos;
			n_actives++;
		}
	}

	Float3 center_of_mass = sum * (1.f / (float)n_actives);
	return center_of_mass;
}
Float3 calcCenterOfMassDistweighted(Row* row, Float3 relative_to) {
	Float3 sum = Float3(0.f);
	int total_weight = 0;

	for (int i = 0; i < row->n_particles; i++) {
		if (!row->particles[i].inactive) {
			float weight = (row->particles[i].pos - relative_to).len();
			sum += row->particles[i].pos * weight;
			total_weight += weight;
		}
	}

	Float3 center_of_mass_distweighted = sum * (1.f / total_weight);
	return center_of_mass_distweighted;
}

void readDataBIN(Row* rows, string path, int N_STEPS, int particles_per_step) {
	char* file_path;
	file_path = &path[0];
	cout << "Reading from file " << file_path << endl;

	FILE* file;
	fopen_s(&file, file_path, "rb");


	uint64_t n_f3s = N_STEPS * particles_per_step * FLOAT3_PER_ATOM;
	printf("Allocating data-buffer of size %.2Lf GB for reading\n", (long double)n_f3s * sizeof(Float3) * 1e-9);
	printf("Here;; %lld\n", n_f3s);
	Float3* buffer = new Float3[n_f3s];
	fread(buffer, sizeof(Float3), n_f3s, file);
	fclose(file);

	printf("Got here");
	for (int step = 0; step < N_STEPS; step++) {
		if (!(step % 100))
			printf("\rReading row %d", step);

		for (int p = 0; p < particles_per_step; p++) {
			uint64_t buffer_index = p * FLOAT3_PER_ATOM + step * particles_per_step * FLOAT3_PER_ATOM;

			rows[step].particles[p] = Particle(p, buffer[buffer_index + 0], buffer[buffer_index + 1]);	// Need of way of determining whether a particle exists, or is in the buffer!!!!!!!!!!!!!!!!!!!!!
		}
	}

	delete[] buffer;
	printf("\n");
}

void markInactiveParticles(Row* row0, Row* row1) {
	for (int i = 0; i < row0->n_particles; i++) {
		if (row0->particles[i].pos.len() == 0 && row1->particles[i].pos.len() == 0) {
			row0->particles[i].inactive = true;
			row1->particles[i].inactive = true;
		}
	}
}

void alignFrameCenters(Row* row0, Row* row1) {
	Float3 com0 = calcCenterOfMass(row0);
	Float3 com1 = calcCenterOfMass(row0);

	Float3 dif = com0 - com1;
	for (int i = 0; i < row0->n_particles; i++) {
		if (!row0->particles[i].inactive) {
			row1->particles[i].pos += dif;
		}
	}		
	// Now both frames share the same com

	for (int i = 0; i < row0->n_particles; i++) {
		if (!row0->particles[i].inactive) {
			row0->particles[i].pos += com0 * -1.f;
			row1->particles[i].pos += com0 * -1.f;
		}
	}
	// Now both frames have a com of (0,0,0)
}


void pitchUp(Row* row) {
	Float3 dir = calcCenterOfMassDistweighted(row, Float3(0.f)).norm();
	Float3 zvector = Float3(0, 0, 1);

	float angle = Float3::getAngle(dir, zvector);
	printf("Pitch angle  %f\n", angle);
	Float3 rot_vector = (zvector.cross(dir)).norm();

	for (int i = 0; i < row->n_particles; i++) {
		if (!row->particles[i].inactive) {
			row->particles[i].pos = Float3::rodriguesRotatation(row->particles[i].pos, rot_vector, -angle);	// RodRot rotates along the FIRST vectors direction, thus we need the minus here
		}
	}
	calcCenterOfMass(row).print('C');
	calcCenterOfMassDistweighted(row, Float3(0,0,0)).print('W');
}

void rotateAroundZ(Row* row, float angle) {
	for (int i = 0; i < row->n_particles; i++) {
		if (!row->particles[i].inactive) {
			row->particles[i].pos = Float3::rodriguesRotatation(row->particles[i].pos, Float3(0,0,1), angle);	
		}
	}
}

void alignPitch(Row* row0, Row* row1) {
	// Assumed that both have a com of 0 at this point!
	pitchUp(row0);
	pitchUp(row1);
	if ((calcCenterOfMassDistweighted(row0, Float3(0.f)).norm() - calcCenterOfMassDistweighted(row1, Float3(0.f)).norm()).len() > 0.1) {
		printf("Error in alignPitch");
		exit(0);
	}		
}

float calcRMSD(Row* row0, Row* row1) {
	double sd = 0;
	int count = 0;
	for (int i = 0; i < row0->n_particles; i++) {
		if (!row0->particles[i].inactive) {
			float l = (row0->particles[i].pos - row1->particles[i].pos).len();
			sd += l*l;		// Dist squared.
			count++;			
		}
	}
	double msd = sd / (double)count;
	double rmsd = sqrt(msd);
	return rmsd;
}

Row* getRotatedCopy(Row* row, float angle) {
	Row* row_ = new Row(row->n_particles);
	for (int i = 0; i < row_->n_particles; i++)
		row_->particles[i] = row->particles[i];
	
	rotateAroundZ(row_, angle);
	return row_;
}

float recurseAlignTilt(Row* row0, Row* row1, float rot_lower, float rot_upper, int depth) {	// Returns optimal angle (lb=lowerbound)
	float rot_middle = (rot_lower + rot_upper) / 2.f;
	//printf("Rot %f %f %f\n",rot_lower, rot_middle, rot_upper);
	float rot_vals[3] = { rot_lower, rot_middle, rot_upper };

	if (depth == 12) {
		return rot_middle;
	}

	float rmsd_vals[3];
	for (int i = 0; i < 3; i++) {		
		Row* rot_frame = getRotatedCopy(row1, rot_vals[i]);
		rmsd_vals[i] = calcRMSD(row0, rot_frame);
		delete rot_frame;
	}


	if (rmsd_vals[0] + rmsd_vals[1] < rmsd_vals[1] + rmsd_vals[2]) {
		return recurseAlignTilt(row0, row1, rot_lower, rot_middle, depth + 1);
	}
	else {
		return recurseAlignTilt(row0, row1, rot_middle, rot_upper, depth + 1);
	}
}


void alignTilt(Row* row0, Row* row1) {
	float rmsd_vals[4];
	float rot_vals[4];
	for (int i = 0; i < 4; i++) {
		float rot = PI * (float)i / 4.f;
		Row* rot_frame = getRotatedCopy(row1, rot);
		rmsd_vals[i] = calcRMSD(row0, rot_frame);
		rot_vals[i] = rot;
		delete rot_frame;
	}
	for (int i = 0; i < 4; i++) {
		for (int ii = i; ii < 4; ii++) {
			if (rmsd_vals[ii] < rmsd_vals[i]) {
				swap(rmsd_vals[i], rmsd_vals[ii]);
				swap(rot_vals[i], rot_vals[ii]);
			}
		}
	}


	float optimal_tilt = recurseAlignTilt(row0, row1, rot_vals[0], rot_vals[1], 0);
	printf("Optimal tilt %f\n\n", optimal_tilt);
	rotateAroundZ(row1, optimal_tilt);
}



void alignFrames(Row* row0, Row* row1) {
	markInactiveParticles(row0, row1);
	alignFrameCenters(row0, row1);
	alignPitch(row0, row1);
	alignTilt(row0, row1);
}

template <typename T>
void printArr(T arr, int n, string arrname) {
	cout << arrname << " = [" << arr[0];
	for (int i = 1; i < n; i++) {
		cout << ',' << arr[i];
	}
	cout << "];\n";
}



int main() {

	const int N_STEPS = 350000;	// Determines file to read
	string workdir = "C:\\PROJECTS\\Quantom\\Simulation\\Steps_" + to_string(N_STEPS) + "\\";
	printf("Hello world");

	int n_compounds = 13;
	int particles_per_compound = 128;
	int particles_per_step = particles_per_compound * n_compounds;

	printf("Allocating %.02Lf GB of RAM for rows\n", (long double)Row::byteSize(particles_per_step) * N_STEPS * 1e-9);	// Gives approx size 3 = 1 in row, 2 in scdp
	Row* rows = new Row[N_STEPS];
	for (int i = 0; i < N_STEPS; i++)
		rows[i] = Row(particles_per_step);

	string path_in = workdir + "sim_traindata.bin";
	readDataBIN(rows, path_in, N_STEPS, particles_per_step);


	const int steps_per_RMSD = 1000;
	const int n_values = N_STEPS / steps_per_RMSD-1;
	int steps[n_values];
	float rmsd[n_values];
	for (int i = 1; i * steps_per_RMSD < N_STEPS; i++) {
		int step = i * steps_per_RMSD;
		alignFrames(&rows[0], &rows[step]);

		steps[i-1] = step;
		rmsd[i-1] = calcRMSD(&rows[0], &rows[step]);
	}


	
	printArr(steps, n_values, "steps");
	printArr(rmsd, n_values, "rmsd");






	return 0;
}