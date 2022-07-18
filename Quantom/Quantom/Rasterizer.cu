#include "Rasterizer.cuh"







void mergeSortAPI(RenderBall* balls, int n_balls);


RenderBall* Rasterizer::render(Simulation* simulation) {    
    solvent_offset = simulation->n_compounds * MAX_COMPOUND_PARTICLES;
    n_threadblocks = (int) ceil((float)simulation->total_particles_upperbound / (float)RAS_THREADS_PER_BLOCK);


	RenderAtom* atoms = getAllAtoms(simulation);


    RenderBall* balls = processAtoms(atoms, simulation);
    mergeSortAPI(balls, simulation->total_particles_upperbound);

	return balls;
}











__global__ void loadCompoundatomsKernel(Box* box, RenderAtom* atoms);
__global__ void loadSolventatomsKernel(Box* box, RenderAtom* atoms, int offset);
__global__ void processAtomsKernel(RenderAtom* atoms, RenderBall* balls);
const int THREADS_PER_LOADSOLVENTSATOMSKERNEL = 100;

RenderAtom* Rasterizer::getAllAtoms(Simulation* simulation) {
	

	RenderAtom* atoms;
	cudaMalloc(&atoms, sizeof(RenderAtom) * simulation->total_particles_upperbound);

    int solvent_blocks = (int) ceil((float)simulation->n_solvents / (float)THREADS_PER_LOADSOLVENTSATOMSKERNEL);


	Box* box = simulation->box;
    if (simulation->n_compounds > 0)
	    loadCompoundatomsKernel << <simulation->n_compounds, MAX_COMPOUND_PARTICLES >> > (box, atoms);
    if (simulation->n_solvents > 0)
    	loadSolventatomsKernel << < solvent_blocks, THREADS_PER_LOADSOLVENTSATOMSKERNEL >> > (simulation->box, atoms, solvent_offset);
	cudaDeviceSynchronize();

	return atoms;
}

void Rasterizer::sortAtoms(RenderAtom* atoms, int dim) {

    cudaDeviceSynchronize();
}

RenderBall* Rasterizer::processAtoms(RenderAtom* atoms, Simulation* simulation) {
    RenderBall* balls_device;
    cudaMalloc(&balls_device, sizeof(RenderBall) * simulation->total_particles_upperbound);
    processAtomsKernel <<< n_threadblocks, RAS_THREADS_PER_BLOCK >>> (atoms, balls_device);
    cudaDeviceSynchronize();

    RenderBall* balls_host = new RenderBall[simulation->total_particles_upperbound];
    cudaMemcpy(balls_host, balls_device, sizeof(RenderBall) * simulation->total_particles_upperbound, cudaMemcpyDeviceToHost);

    cudaFree(atoms);
    cudaFree(balls_device);

    return balls_host;
}







__device__ ATOM_TYPE RAS_getTypeFromIndex(int atom_index) {
    switch (atom_index)
    {
    case 0:
        return ATOM_TYPE::SOL;
    case 1:
        return ATOM_TYPE::C;
    case 2:
        return ATOM_TYPE::O;
    case 3:
        return ATOM_TYPE::N;
    case 4:
        return ATOM_TYPE::H;
    case 5: 
        return ATOM_TYPE::P;
    case 6:
        return ATOM_TYPE::SOL;
    default:
        return ATOM_TYPE::NONE;
    }
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









__global__ void loadCompoundatomsKernel(Box * box, RenderAtom * atoms) {                                                            // TODO: CAN ONLY HANDLE ONE COMPOUND!!!
    int local_id = threadIdx.x;
    int compound_id = blockIdx.x;
    int global_id = threadIdx.x + blockIdx.x * blockDim.x;

    
    if (local_id < box->compounds[compound_id].n_particles) {
        atoms[global_id].pos = box->compound_state_array[compound_id].positions[local_id];                                                          // Might need to change this, if thread> n_particles!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //atoms[global_id].pos.print('A');
        atoms[global_id].mass = SOLVENT_MASS;                                                         // TEMP
        //atoms[global_id].atom_type = RAS_getTypeFromIndex(box->compounds[compound_id].atom_types[local_id]);
        atoms[global_id].atom_type = RAS_getTypeFromIndex(box->compounds[compound_id].atom_color_types[local_id]);
    }
    else {
        atoms[global_id].atom_type = ATOM_TYPE::NONE;
    }
}

__global__ void loadSolventatomsKernel(Box* box, RenderAtom * atoms, int offset) {
    int solvent_id = threadIdx.x + blockIdx.x * THREADS_PER_LOADSOLVENTSATOMSKERNEL;

    if (solvent_id < box->n_solvents) {
        atoms[solvent_id + offset].pos = box->solvents[solvent_id].pos;
        atoms[solvent_id + offset].mass = SOLVENT_MASS;
        atoms[solvent_id + offset].atom_type = SOL;
    }    
}




__global__ void processAtomsKernel(RenderAtom* atoms, RenderBall* balls) { 
    int index = threadIdx.x + blockIdx.x * RAS_THREADS_PER_BLOCK;
    

    RenderAtom atom = atoms[index];

    atom.color = getColor(atom.atom_type);
    atom.radius = (getRadius(atom.atom_type)) / (1.f+atom.pos.y * 0.001f);       // [nm]

    //atoms[index] = atom;

    // Convert units to normalized units for OpenGL
    atom.radius = 0.25 * atom.radius;            // Yeah, i'm just eyeballing this..
    for (int dim = 0; dim < 3; dim++) {
        *atom.pos.placeAt(dim) = (atom.pos.at(dim) / (double) BOX_LEN - 0.5f) * 1.8l;
    }


    RenderBall ball(atom.pos, atom.radius, atom.color);
    if (atom.atom_type == ATOM_TYPE::NONE)
        ball.disable = true;
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


    RenderBall* sorted_balls = mergeSort(balls, 0, n_balls - 1);
    for (int i = 0; i < n_balls; i++) {
        balls[i] = sorted_balls[i];
    }

    int* mapping = new int[n_balls];
    //for (int i = 0; i < n_atoms; i++) mapping[i] = sorted_atoms[i].id;

    delete[] sorted_balls;
}