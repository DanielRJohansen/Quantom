#include "Rasterizer.cuh"







void mergeSortAPI(RenderBall* balls, int n_balls);


RenderBall* Rasterizer::render(Simulation* simulation) {    
    actual_n_particles = simulation->box->compounds[0].n_particles + simulation->n_solvents;
    solvent_offset = simulation->box->compounds[0].n_particles;
    n_threadblocks = ceil((float)actual_n_particles / (float)RAS_THREADS_PER_BLOCK);


	RenderAtom* atoms = getAllAtoms(simulation);

    //sortAtoms(atoms, 1);    // dim 1 is y

    RenderBall* balls = processAtoms(atoms);
    mergeSortAPI(balls, actual_n_particles);

	return balls;
}











__global__ void loadCompoundatomsKernel(Box* box, RenderAtom* atoms);
__global__ void loadSolventatomsKernel(Box* box, RenderAtom* atoms, int offset);
__global__ void processAtomsKernel(RenderAtom* atoms, RenderBall* balls, int n_atoms);


RenderAtom* Rasterizer::getAllAtoms(Simulation* simulation) {
	

	RenderAtom* atoms;
	cudaMalloc(&atoms, sizeof(RenderAtom) * actual_n_particles);

	Box* box = simulation->box;
	loadCompoundatomsKernel << <1, simulation->box->compounds[0].n_particles >> > (box, atoms);
	loadSolventatomsKernel << <1, simulation->n_solvents >> > (simulation->box, atoms, solvent_offset);
	cudaDeviceSynchronize();

	return atoms;
}

void Rasterizer::sortAtoms(RenderAtom* atoms, int dim) {

    cudaDeviceSynchronize();
}

RenderBall* Rasterizer::processAtoms(RenderAtom* atoms) {
    RenderBall* balls_device;
    cudaMalloc(&balls_device, sizeof(RenderBall) * actual_n_particles);

    processAtomsKernel <<< n_threadblocks, RAS_THREADS_PER_BLOCK >>> (atoms, balls_device, actual_n_particles);
    cudaDeviceSynchronize();

    RenderBall* balls_host = new RenderBall[actual_n_particles];
    cudaMemcpy(balls_host, balls_device, sizeof(RenderBall) * actual_n_particles, cudaMemcpyDeviceToHost);

    cudaFree(atoms);
    cudaFree(balls_device);

    return balls_host;
}













__device__ ATOM_TYPE RAS_getTypeFromMass(double mass) {
    mass *= 1000.f;   //convert to g
    if (mass < 4)
        return ATOM_TYPE::H;
    if (mass < 14)
        return ATOM_TYPE::C;
    if (mass < 15)
        return ATOM_TYPE::N;
    if (mass < 18)
        return ATOM_TYPE::O;
    if (mass < 32)
        return ATOM_TYPE::P;
    return ATOM_TYPE::NONE;
}

__device__ Int3 getColor(ATOM_TYPE atom_type) {
    switch (atom_type)
    {
    case ATOM_TYPE::SOL:
        return Int3(0x03, 0xa9, 0xf4);
    case ATOM_TYPE::H:
        return Int3(0xF1, 0xF1, 0xF1);
    case ATOM_TYPE::O:
        return Int3(0xE0, 0x20, 0x20);
    case ATOM_TYPE::C:
        return Int3(0x30, 0x10, 0x90);
    case ATOM_TYPE::P:
        return Int3(0xFC, 0xF7, 0x5E);
    case ATOM_TYPE::N:
        return Int3(0x2E, 0x8B, 0x57);
    case ATOM_TYPE::NONE:
        return Int3(0xF2, 0xE5, 0xD9);     
    default:
        return Int3(0, 0, 0);
    }
}

__device__ float getRadius(ATOM_TYPE atom_type) {
    switch (atom_type)
    {

    case ATOM_TYPE::H:
        return 0.04;
    case ATOM_TYPE::C:
        return 0.1;
    case ATOM_TYPE::N:
        return 0.1;
    case ATOM_TYPE::O:
        return 0.12;
    case ATOM_TYPE::SOL:
        return 0.08;
    case ATOM_TYPE::P:
        return 0.15;
    case ATOM_TYPE::NONE:
        return 1;
    default:
        return 1;
    }
}









__global__ void loadCompoundatomsKernel(Box * box, RenderAtom * atoms) {
    atoms[threadIdx.x].pos = box->compound_state_array[0].positions[threadIdx.x];
    atoms[threadIdx.x].mass = box->compounds[0].particles[threadIdx.x].mass;
}

__global__ void loadSolventatomsKernel(Box * box, RenderAtom * atoms, int offset) {
    atoms[threadIdx.x + offset].pos = box->solvents[threadIdx.x].pos;
    atoms[threadIdx.x + offset].mass = SOLVENT_MASS;
    atoms[threadIdx.x + offset].atom_type = SOL;
}

__global__ void processAtomsKernel(RenderAtom* atoms, RenderBall* balls, int n_atoms) {
    int index = threadIdx.x + blockIdx.x * RAS_THREADS_PER_BLOCK;
    
    if (index > n_atoms)
        return;

    RenderAtom atom = atoms[index];

    atom.atom_type = atom.atom_type == NONE ? RAS_getTypeFromMass(atom.mass) : atom.atom_type;

    atom.color = getColor(atom.atom_type);
    atom.radius = (getRadius(atom.atom_type)) / (1.f+atom.pos.y * 0.001);       // [nm]

    //atoms[index] = atom;

    // Convert units to normalized units for OpenGL
    atom.radius = 0.5 * atom.radius;            // Yeah, i'm just eyeballing this..
    for (int dim = 0; dim < 3; dim++) {
        *atom.pos.placeAt(dim) = (atom.pos.at(dim) / BOX_LEN - 0.5f) * 1.8;
    }


    RenderBall ball(atom.pos, atom.radius, atom.color);
    balls[index] = ball;
}


RenderBall* merge(const RenderBall* left, int n_left, const RenderBall* right, int n_right) {
    RenderBall* merged = new RenderBall[n_left + n_right];
    int l = 0;
    int r = 0;
    int index = 0;

    while (l < n_left && r < n_right) {
        if (left[l].pos.y < right[r].pos.y) {
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

RenderBall* mergeSort(const RenderBall* atoms, const int l, const int r) {	// l and r are indexes of the two extremes
    int m = (l + r) / 2;

    RenderBall* left;
    RenderBall* right;

    //printf("step\n");
    if (r - l > 1) {
        //printf("l %d r %d\n", l, r);
        left = mergeSort(atoms, l, m - 1);
        right = mergeSort(atoms, m, r);
        RenderBall* merged = merge(left, m - l, right, r - m + 1);

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
        return new RenderBall(atoms[l]);
    }
}

void mergeSortAPI(RenderBall* balls, int n_balls) {					// Returns a mapping where mapping[0] is the closest id, mapping [1] seconds closest, so on

    //RenderBall* balls_temp = new RenderBall[n_balls];
    //memcpy(balls_temp, balls, sizeof(RenderBall) * n_balls);




    RenderBall* sorted_balls = mergeSort(balls, 0, n_balls - 1);
    for (int i = 0; i < n_balls; i++) {
        balls[i] = sorted_balls[i];
    }

    int* mapping = new int[n_balls];
    //for (int i = 0; i < n_atoms; i++) mapping[i] = sorted_atoms[i].id;

    delete[] sorted_balls;
}