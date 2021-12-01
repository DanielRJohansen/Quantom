#include "Raytracer.cuh"



Ray::Ray(Float3 unit_vector, Float3 origin) : unit_vector(unit_vector), origin(origin) {}

void Ray::reset() {
    closest_collision = 999999999;
    atom_type = -1;
    log_particle = false;
}

__device__ bool Ray::hitsParticle(Float3* particle_center, double particle_radius) {
    return (distToPoint(*particle_center) < particle_radius);
}

__device__ double Ray::distToSphereIntersect(Float3* particle_center, double particle_radius) {
    Float3 projection_on_ray = origin + unit_vector * ((*particle_center - origin).dot(unit_vector) / unit_vector.dot(unit_vector));
    double center_to_projection = (projection_on_ray - *particle_center).len();
    double projection_to_intersect = sqrtf(particle_radius * particle_radius - center_to_projection * center_to_projection);
    return (projection_on_ray - origin).len() - projection_to_intersect;
}


__device__ double Ray::distToPoint(Float3 point) {
    Float3 far_ray_point = origin + unit_vector * 99999999;	////BAAAAAAAAAAAAAD
    return (
        ((far_ray_point - origin).cross(origin - point)).len()
        /
        (far_ray_point - origin).len()
        );
}
    
__device__ void Ray::searchCompound(CompoundState* state, Box* box, int compound_index) {
    for (int i = 0; i < state->n_particles; i++) {
        if (hitsParticle(&state->positions[i], 0.170)) {              // LOOK HERE at index 0....
            double dist = distToSphereIntersect(&state->positions[i], 0.170);
            if (dist < closest_collision) {
                closest_collision = dist;
                atom_type = 0;
                //if (compound_index == LOGBLOCK && i == LOGTHREAD)
                  //  log_particle = true;        // Temp
            }
        }
    }
}

__device__ bool Ray::searchParticle(Float3* pos, int index, bool is_virtual)
{
    if (hitsParticle(pos, 0.150)) {
        double dist = distToSphereIntersect(pos, 0.150);
        if (dist < closest_collision) {
            closest_collision = dist;
            atom_type = 1;
            if ((index == LOGTHREAD && LOGTYPE == 0) || is_virtual)
                log_particle = true;
            return true;
        }
    }
    return false;
}


__device__ void paintImage(uint8_t* image, Ray* ray) {
    if (ray->atom_type == 0) {      // Carbon
        image[0] = 40;
        image[1] = 0;
        image[2] = 40;
        image[3] = 255;
    }
    else if (ray->atom_type == 1) {      // Solvent
        image[0] = 200;
        image[3] = 255;
    }
    else {
        image[0] = 0xFE;
        image[1] = 0xFE;
        image[2] = 0xFA;
        image[3] = 0xE2;
    }

    if (ray->log_particle) {
        image[0] = 0x53;
        image[1] = 0xAF;
        image[2] = 0x8B;
    }
}



Raytracer::Raytracer(Simulation* simulation, bool verbose) {
    auto t0 = std::chrono::high_resolution_clock::now();
    printf("\n\n");
    setGPU();


    double base = 0;
	double principal_point_increment = (double) (BOX_LEN) / (double)RAYS_PER_DIM;

	Ray* host_rayptr = new Ray[NUM_RAYS];
	focalpoint = Float3(BOX_LEN/2.f, -(BOX_LEN / 2.f) * FOCAL_LEN_RATIO - BOX_LEN, BOX_LEN / 2.f);

    

	int index = 0;
    for (int z_index = 0; z_index < RAYS_PER_DIM; z_index++) {
        for (int x_index = 0; x_index < RAYS_PER_DIM; x_index++) {
            double z = base + principal_point_increment * (double)z_index;
            double x = base + principal_point_increment * (double)x_index;
			Float3 vector = Float3(x, base, z) - focalpoint;
            host_rayptr[index++] = Ray(vector.norm(), focalpoint);
		}
	}

    cuda_status = cudaMallocManaged(&rayptr, NUM_RAYS * sizeof(Ray));
    cuda_status = cudaMemcpy(rayptr, host_rayptr, NUM_RAYS * sizeof(Ray), cudaMemcpyHostToDevice);

    if (verbose) {
        printf("Focal point: %.3f %.3f %.3f\n", focalpoint.x, focalpoint.y, focalpoint.z);
        printf("Allocating %d MB of ram for Rays... \n", NUM_RAYS * sizeof(Ray) / 1000000);
    }
    

    int duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t0).count();
    printf("\nRaytracer initiated in %d ms\n\n", duration);
}





__global__ void renderKernel(Ray* rayptr, uint8_t* image, Box* box) {
    int  index = blockIdx.x * blockDim.x + threadIdx.x;
    Ray ray = rayptr[index];
    ray.reset();
    
    

    for (int i = 0; i < box->n_compounds; i++) {
        CompoundState* compoundstate = &box->compound_state_array[i];
        ray.searchCompound(compoundstate, box, i);
    }

    for (int i = 0; i < box->n_solvents; i++) {
        if (ray.searchParticle(&box->solvents[i].pos, i)) {
            
        }
    }
    
    
    paintImage(&image[index * 4], &ray);
}

__global__ void renderKernel(Ray* rayptr, uint8_t* image, Trajectory* trajectory, int step) {
    int  index = blockIdx.x * blockDim.x + threadIdx.x;
    Ray ray = rayptr[index];
    ray.reset();
    
    int step_offset = trajectory->n_particles * step;
    for (int i = 0; i < trajectory->n_particles; i++) {
        Float3* particle_pos = &trajectory->positions[step_offset + i];
        bool is_virtual = (i == 1);
        if (ray.searchParticle(particle_pos, i, is_virtual)) {

        }
    }

    paintImage(&image[index * 4], &ray);
}







    
uint8_t* Raytracer::render(Simulation* simulation) {
    auto start = std::chrono::high_resolution_clock::now();

    cuda_status = cudaGetLastError();
    if (cuda_status != cudaSuccess) {
        fprintf(stderr, "Something is wrong");
        exit(1);
    }



    cudaStream_t renderstream;
    cudaStreamCreate(&renderstream);


    uint8_t* cuda_image;
    int im_bytesize = NUM_RAYS * 4 * sizeof(uint8_t);
    cudaMallocManaged(&cuda_image, im_bytesize);


    renderKernel << < RAYS_PER_DIM, RAYS_PER_DIM, 0>>> ( rayptr, cuda_image, simulation->box);
    uint8_t* image = new uint8_t[NUM_RAYS * 4];
    cudaMemcpy(image, cuda_image, im_bytesize, cudaMemcpyDeviceToHost);


    cudaFree(cuda_image);
    cudaStreamDestroy(renderstream);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    printf("\tRender time: %4d ms  ", duration.count());
    // First render: 666 ms

    return image;
}

uint8_t* Raytracer::render(Trajectory* trajectory, int step)
{
    uint8_t* cuda_image;
    int im_bytesize = NUM_RAYS * 4 * sizeof(uint8_t);
    cudaMallocManaged(&cuda_image, im_bytesize);


    renderKernel << < RAYS_PER_DIM, RAYS_PER_DIM>> > (rayptr, cuda_image, trajectory, step);
    uint8_t* image = new uint8_t[NUM_RAYS * 4];
    cudaMemcpy(image, cuda_image, im_bytesize, cudaMemcpyDeviceToHost);


    cudaFree(cuda_image);
    return image;
}
