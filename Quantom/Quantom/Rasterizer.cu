#include "Rasterizer.cuh"









RenderBall* Rasterizer::render(Simulation* simulation) {    
    actual_n_particles = simulation->box->compounds[0].n_particles + simulation->n_solvents;
    solvent_offset = simulation->box->compounds[0].n_particles;
    n_threadblocks = ceil((float)actual_n_particles / (float)RAS_THREADS_PER_BLOCK);


	RenderAtom* atoms = getAllAtoms(simulation);

    sortAtoms(atoms, 1);    // dim 1 is y

    RenderBall* balls = processAtoms(atoms);


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
