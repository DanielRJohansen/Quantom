#include "Raytracer.cuh"



Ray::Ray(Float3 unit_vector, Float3 origin) : unit_vector(unit_vector), origin(origin) {}

void Ray::reset() {
    closest_collision = 999999999;
    illumination = 1.f;
    atom_type = NONE;
    log_particle = false;
}


__device__ ATOM_TYPE getTypeFromMass(double mass) {
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
    
__device__ int Ray::searchCompound(CompoundState* state, Box* box, int compound_index) {    // Returns compound-local index of hit particle
    int index = -1;
    for (int i = 0; i < state->n_particles; i++) {
        if (hitsParticle(&state->positions[i], 0.170)) {              // LOOK HERE at index 0....
            double dist = distToSphereIntersect(&state->positions[i], 0.170);
            if (dist < closest_collision) {
                closest_collision = dist;
                index = i;
                //atom_type = ATOM_TYPE::C;
                illumination = (((origin + unit_vector * dist).z - state->positions[i].z) + 0.150) / 0.300;

                //if (compound_index == LOGBLOCK && i == LOGTHREAD)
                  //  log_particle = true;        // Temp
            }
        }
    }
    return index;
}

__device__ bool Ray::searchParticle(Float3* pos, int index, bool is_virtual)
{
    if (hitsParticle(pos, 0.150)) {
        double dist = distToSphereIntersect(pos, 0.150);
        if (dist < closest_collision) {
            closest_collision = dist;
            
            //illumination = (((origin + unit_vector * dist).z - pos->z) + 0.150) / 0.300;   // + 
            illumination = (((origin + unit_vector * dist).z - pos->z) + 0.200) / 0.350;   // Offset by 50 to avoid black atoms
            return true;
        }
    }
    return false;
}


__device__ void paintImage(uint8_t* image, Ray* ray) {
    switch (ray->atom_type)
    {
    case ATOM_TYPE::SOL:
        image[0] = 0x03;
        image[1] = 0xa9;
        image[2] = 0xf4;
        image[3] = 255;
        break;
    case ATOM_TYPE::H:
        image[0] = 250;
        image[1] = 250;
        image[2] = 250;
        image[3] = 255;
        break;
    case ATOM_TYPE::O:
        image[0] = 240;
        image[1] = 20;
        image[2] = 20;
        image[3] = 255;
        break;
    case ATOM_TYPE::C:
        image[0] = 40;
        image[1] = 10;
        image[2] = 180;
        image[3] = 255;
        break;
    case ATOM_TYPE::P:
        image[0] = 0xFC;
        image[1] = 0xF7;
        image[2] = 0x5E;
        image[3] = 0xFF;
        break;
    case ATOM_TYPE::N:
        image[0] = 0x2E;
        image[1] = 0x8B;
        image[2] = 0x57;
        image[3] = 0xFF;
        break;
    case ATOM_TYPE::NONE:
        image[0] = 0xF2;
        image[1] = 0xE5;
        image[2] = 0xD9;
        image[3] = 0xFF;
        break;
    default:
        image[3] = 0x00;
        break;
    }

    if (ray->illumination > 1.f)
        ray->illumination = 1;
    image[3] *= ray->illumination;

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
    
    

    for (int compound_index = 0; compound_index < box->n_compounds; compound_index++) {
        CompoundState* compoundstate = &box->compound_state_array[compound_index];
        int particle_index = ray.searchCompound(compoundstate, box, compound_index);
        if (particle_index != -1) {            
            ray.atom_type = getTypeFromMass(box->compounds[compound_index].particles[particle_index].mass);
        }
    }

    for (int i = 0; i < box->n_solvents; i++) {
        if (ray.searchParticle(&box->solvents[i].pos, i)) {
            ray.atom_type = ATOM_TYPE::SOL;
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



    uint8_t* cuda_image;
    int im_bytesize = NUM_RAYS * 4 * sizeof(uint8_t);
    cudaMallocManaged(&cuda_image, im_bytesize);

    cudaDeviceSynchronize();
    renderKernel << < RAYS_PER_DIM, RAYS_PER_DIM, 0>>> ( rayptr, cuda_image, simulation->box);
    cudaDeviceSynchronize();

    uint8_t* image = new uint8_t[NUM_RAYS * 4];
    cudaMemcpy(image, cuda_image, im_bytesize, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();

    cudaFree(cuda_image);
    //cudaStreamDestroy(renderstream);
    cudaDeviceSynchronize();

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    printf("\tRender time: %4d ms  ", duration.count());

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
