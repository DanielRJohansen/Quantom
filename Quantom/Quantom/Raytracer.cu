#include "Raytracer.cuh"



Ray::Ray(Float3 unit_vector, Float3 origin) : unit_vector(unit_vector), origin(origin) {}

void Ray::reset() {
    closest_collision = 999999999;
    atom_type = -1;
}

__device__ bool Ray::hitsParticle(Float3* particle_center, float particle_radius) {
    //printf("dist: %f \n", distToPoint(*particle_center));
    return (distToPoint(*particle_center) < particle_radius);
}

__device__ float Ray::distToSphereIntersect(Float3* particle_center, float particle_radius) {
    Float3 projection_on_ray = origin + unit_vector * ((*particle_center - origin).dot(unit_vector) / unit_vector.dot(unit_vector));
    float center_to_projection = (projection_on_ray - *particle_center).len();
    float projection_to_intersect = sqrtf(particle_radius * particle_radius - center_to_projection * center_to_projection);
    return (projection_on_ray - origin).len() - projection_to_intersect;
}


__device__ float Ray::distToPoint(Float3 point) {
    Float3 far_ray_point = origin + unit_vector * 99999999;	////BAAAAAAAAAAAAAD
    return (
        ((far_ray_point - origin).cross(origin - point)).len()
        /
        (far_ray_point - origin).len()
        );
}
    
__device__ void Ray::searchCompound(CompoundState* state, Box* box) {
    for (int i = 0; i < state->particle_cnt; i++) {
        if (hitsParticle(&state->positions[i], box->rendermolecule.radii[i])) {
            if (distToSphereIntersect(&state->positions[i], box->rendermolecule.radii[i]) < closest_collision) {
                atom_type = i;
            }
        }
    }
}







Raytracer::Raytracer(Simulation* simulation, bool verbose) {
    auto t0 = std::chrono::high_resolution_clock::now();
    printf("\n\n");
    setGPU();


    float base = 0;
	float principal_point_increment = (float) (BOX_LEN) / (float)RAYS_PER_DIM;

	Ray* host_rayptr = new Ray[NUM_RAYS];
	focalpoint = Float3(BOX_LEN/2.f, -(BOX_LEN / 2.f) * FOCAL_LEN_RATIO - BOX_LEN, BOX_LEN / 2.f);

    

	int index = 0;
    for (int z_index = 0; z_index < RAYS_PER_DIM; z_index++) {
        for (int x_index = 0; x_index < RAYS_PER_DIM; x_index++) {
            float z = base + principal_point_increment * (float)z_index;
            float x = base + principal_point_increment * (float)x_index;
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
        CompoundState* compoundstate = &box->compound_state_buffer[i];
        ray.searchCompound(compoundstate, box);
    }


    for (int i = 0; i < 3; i++) {
        image[index * 4 + i] = box->rendermolecule.colors[ray.atom_type][i];
    }
    if (ray.atom_type != -1)
        image[index * 4 + 3] = 255;
    else
        image[index * 4 + 3] = 0;
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