#include "Raytracer.cuh"



Ray::Ray(Float3 unit_vector, Float3 origin) : unit_vector(unit_vector), origin(origin) {}

__device__ void Ray::findBlockHits(Box* box, Float3 focalpoint) {
    for (int i = 0; i < MAX_RAY_BLOCKS; i++)
		block_indexes[i] = -1;

    int block_index = 0;
    int bpd = box->blocks_per_dim;

    for (int y = 0; y < bpd; y++) {
        for (int z = 0; z < bpd; z++) {
            for (int x = 0; x < bpd; x++) {
                if (block_index == MAX_RAY_BLOCKS)
                    break;

                int index = z * bpd * bpd + y * bpd + x;


                if (hitsBlock(&box->blocks[index], focalpoint)) {
                    block_indexes[block_index] = index;

                    block_index++;
                    
                        
                }
            }
        }
    }

    
	for (int i = 1; i < block_index; i++) {		// I think only 1 swap will be necessary...
		float d1_sq = (box->blocks[block_indexes[i - 1]].center - focalpoint).lenSquared();
		float d2_sq = (box->blocks[block_indexes[i]].center - focalpoint).lenSquared();

		if (d2_sq < d1_sq) {
			int tmp = block_indexes[i - 1];
			block_indexes[i - 1] = block_indexes[i];
			block_indexes[i] = tmp;
		}
	}
}

__device__ float cudaMax(float a, float b) { 
    if (a > b)
        return a;
    return b;
}
__device__ float cudaMin(float a, float b) {
    if (a < b)
        return a;
    return b;
}

__device__ bool containedBetweenPlanes(float plane_min, float plane_max, float dir_dim, float origin_dim) {
    float tmin = -999999;
    float tmax = 999999;

    float invD = 1.0f / dir_dim;
    float t0 = (plane_min - origin_dim) * invD;
    float t1 = (plane_max - origin_dim) * invD;

    if (invD < 0.0f) {
        float temp = t1;
        t1 = t0;
        t0 = temp;
    }

    tmin = t0 > tmin ? t0 : tmin;
    tmax = t1 < tmax ? t1 : tmax;



    if (tmax < tmin)    // was <=
        return false;

    return true;
}

__device__ bool Ray::hitsBlock(Block* block, Float3 focalpoint) {
    float a = BLOCK_LEN;               // FUCK YOU VISUAL STUDIO
    Float3 blocksize = Float3(a, a, a);

    Float3 offset(BODY_RADIUS, BODY_RADIUS, BODY_RADIUS);   // Radius of largest molecule!
    Float3 min = block->center - Float3(FOCUS_LEN_HALF, FOCUS_LEN_HALF, FOCUS_LEN_HALF) - offset *2;    //I have no idea why we need *2
    Float3 max = block->center + Float3(FOCUS_LEN_HALF, FOCUS_LEN_HALF, FOCUS_LEN_HALF) + offset *2;

    float tmin = -DBL_MAX;
    float tmax = DBL_MAX;

    for (int dim = 0; dim < 3; dim++) {
        float invD = unit_vector.at(dim) != 0 ? 1.0f / unit_vector.at(dim) : 999999999;

        float t0 = (min.at(dim) - focalpoint.at(dim)) * invD;
        float t1 = (max.at(dim) - focalpoint.at(dim)) * invD;

        if (invD < 0.0f) {
            float temp = t1;
            t1 = t0;
            t0 = temp;
        }

        tmin = t0 > tmin ? t0 : tmin;
        tmax = t1 < tmax ? t1 : tmax;


        if (tmax <= tmin) {    // was <=
            return false;
        }
    } 

        
    return true;
}

__device__ bool Ray::hitsBlock(Float3* blockmin, Float3* blockmax, Float3* focalpoint) {
    float tmin = -DBL_MAX;
    float tmax = DBL_MAX;
    
    for (int dim = 0; dim < 3; dim++) {
        float invD = unit_vector.at(dim) != 0 ? 1.0f / unit_vector.at(dim) : 999999999;

        float t0 = (blockmin->at(dim) - focalpoint->at(dim)) * invD;
        float t1 = (blockmax->at(dim) - focalpoint->at(dim)) * invD;

        if (invD < 0.0f) {
            float temp = t1;
            t1 = t0;
            t0 = temp;
        }

        tmin = t0 > tmin ? t0 : tmin;
        tmax = t1 < tmax ? t1 : tmax;


        if (tmax <= tmin) {    // was <=
            return false;
        }
    }


    return true;
}

__device__ bool Ray::hitsParticle(Particle* particle) {
    //if (distToPoint(particle->pos) < BODY_RADIUS)
    if (distToPoint(particle->pos) < particle->radius)
        return true;
    return false;
}

__device__ bool Ray::searchBlock(Block* block, uint8_t* image) {
    Particle* closest_particle;
    float big_number = 99999;
    float closest_intersect = big_number;
    int  index = blockIdx.x * blockDim.x + threadIdx.x;


    for (int j = 0; j < MAX_FOCUS_BODIES; j++) {
        Particle* particle = &block->focus_particles[j];
        if (!particle->active) {  // We can do this because bodies are loaded from index 0 at each timestep. MIGHT NEED TO CHANGE IN THE FUTURE!
            break;
        }

        if (hitsParticle(&block->focus_particles[j])) {
            
            float dist = distToSphereIntersect(particle);
            if (dist < closest_intersect) {
                closest_intersect = dist;
                closest_particle = particle;
            }
            //if (ray.moleculeCollisionHandling(&block->focus_particles[j], mol_library, image)) {
             //   return;
            //}
        }
    }
    if (closest_intersect < big_number) {
        image[index * 4 + 0] = closest_particle->color[0];
        image[index * 4 + 1] = closest_particle->color[1];
        image[index * 4 + 2] = closest_particle->color[2];
        image[index * 4 + 3] = 255;
        return true;
    }
    return false;
}



__device__ bool Ray::moleculeCollisionHandling(Particle* particle, MoleculeLibrary* mol_library, uint8_t* image) {
    Molecule* mol = &mol_library->molecules[0];

    const int infinity = 9999999;

    float closest_collision = infinity;
    Atom closest_atom;
    closest_atom.pos = Float3(0, infinity, 0);  // Make sure its infinitely far away in y direction.




    //Float3 molecule_tilt_vector = Float3(0, 1, 0).rotateAroundOrigin(body->rotation);
    /*
    for (int atom_index = 0; atom_index < mol->n_atoms; atom_index++) {
        
        // Local copy which we can manipulate
        Atom atom = mol->atoms[atom_index];
        atom.pos = atom.pos.rotateAroundVector(body->rotation, molecule_tilt_vector);   // Rotate around its relative origin, before moving pos to global coords
        atom.pos = atom.pos + body->pos;

        Float3 atom_absolute_pos = body->pos + atom.pos;


        if (distToPoint(atom.pos) < atom.radius) {
            float collision_dist = distToSphereIntersect(&atom);
            if (collision_dist < closest_collision) {
                closest_atom = atom;
                closest_collision = collision_dist;
            }


        }
            
    }
    */

    int index = blockIdx.x * blockDim.x + threadIdx.x;
    

    image[index * 4 + 0] = closest_atom.color[0];
    image[index * 4 + 1] = closest_atom.color[1];
    image[index * 4 + 2] = closest_atom.color[2];
    image[index * 4 + 3] = 255;


    /*  // TEMP REMOVED 
    if (closest_atom.pos.y != infinity) {
        int index = blockIdx.x * blockDim.x + threadIdx.x;

        image[index * 4 + 0] = closest_atom.color[0];
        image[index * 4 + 1] = closest_atom.color[1];
        image[index * 4 + 2] = closest_atom.color[2];
        image[index * 4 + 3] = 255;
        return true;
    }*/
    
    
    return false;
}

/*__device__ float Ray::distToSphereIntersect(Atom* atom) {
    Float3 projection_on_ray = origin + unit_vector * ((atom->pos - origin).dot(unit_vector) / unit_vector.dot(unit_vector));
    float center_to_projection = (projection_on_ray - atom->pos).len();
    float projection_to_intersect = sqrtf(atom->radius * atom->radius - center_to_projection * center_to_projection);
    return (projection_on_ray - origin).len() - projection_to_intersect;
}*/
__device__ float Ray::distToSphereIntersect(Particle* particle) {
    Float3 projection_on_ray = origin + unit_vector * ((particle->pos - origin).dot(unit_vector) / unit_vector.dot(unit_vector));
    float center_to_projection = (projection_on_ray - particle->pos).len();
    float projection_to_intersect = sqrtf(particle->radius * particle->radius - center_to_projection * center_to_projection);
    return (projection_on_ray - origin).len() - projection_to_intersect;
}

__global__ void initRayKernel(Ray* rayptr, Box* box, Float3 focalp, int offset) {
	int index = blockIdx.x * blockDim.x + threadIdx.x + offset;
    //int index = blockIdx.x * stride + threadIdx.x;
    if (index > NUM_RAYS)
        return;

    
    //extern __shared__ Float3 row_blockcenters[];
    extern __shared__ Float3 row_block_minsmaxs[]; //[min0, max0,min1,max1]
    __shared__ int bpd;
    __shared__ Float3 boxsize_half;
    __shared__ Float3 boxsize_offset;
    __shared__ Float3 focalpoint;


    if (threadIdx.x == 0) {
        bpd = box->blocks_per_dim;
        boxsize_half = Float3(FOCUS_LEN_HALF, FOCUS_LEN_HALF, FOCUS_LEN_HALF);
        boxsize_offset = Float3(BODY_RADIUS, BODY_RADIUS, BODY_RADIUS);
        focalpoint = focalp;
    }
    __syncthreads();




    Ray ray = rayptr[index];
    int ray_blockcnt = 0;
    for (int i = 0; i < MAX_RAY_BLOCKS; i++)
        ray.block_indexes[i] = -1;





    for (int z = 0; z < bpd; z++) {
        for (int y = 0; y < bpd; y++) {

            // First load the entire row into shared memory
            if (threadIdx.x < bpd) {            // UNSAFE if i change bpd too large, og threads_per_gpublock too small!!!
                int blockindex1 = z * bpd * bpd + y * bpd + threadIdx.x;
                Float3 block_center = box->blocks[blockindex1].center;
                Float3 min = block_center - boxsize_half - boxsize_offset;// -offset;    //I have no idea why we need *2
                Float3 max = block_center + boxsize_half + boxsize_offset;

                row_block_minsmaxs[threadIdx.x * 2] = min;
                row_block_minsmaxs[threadIdx.x * 2 + 1] = max;
                //row_block_minsmaxs[threadIdx.x] = box->blocks[blockindex1].center;
            }
            __syncthreads();

            
            // Each threads loops over all blockcenters
            for (int x = 0; x < bpd; x++) {
                if (ray_blockcnt == MAX_RAY_BLOCKS)
                    break;  // Only break this inner loop, as the thread might need to fetch blockcenters in other rows

                //break;
                //if (ray.hitsBlock(row_block_minsmaxs[x], focalpoint)) {
                if (ray.hitsBlock(&row_block_minsmaxs[x*2], &row_block_minsmaxs[x * 2 + 1], &focalpoint)) {
                    ray.block_indexes[ray_blockcnt++] = z * bpd * bpd + y * bpd + x;
                }
            }
            


        }
    }
    



	//ray.findBlockHits(box, focalpoint);

    rayptr[index] = ray;
}


    

        
        
        










Raytracer::Raytracer(Simulation* simulation, bool verbose) {
    auto t0 = std::chrono::high_resolution_clock::now();
    printf("\n\n");
    setGPU();


    float base = -(BOX_LEN) / 2.f;
	float principal_point_increment = (BOX_LEN) / (float)RAYS_PER_DIM;

	Ray* host_rayptr = new Ray[NUM_RAYS];
	focalpoint = Float3(0, -(BOX_LEN / 2.f) * FOCAL_LEN_RATIO - BOX_LEN, 0);

    

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
    

    // INIT KERNEL PREP ----------------------------------------------------
    const int n_streams = 10;
    cudaStream_t stream[n_streams];
    for (int i = 0; i < n_streams; i++)
        cudaStreamCreate(&stream[i]);
    cudaDeviceSynchronize();
    cuda_status = cudaGetLastError();
    if (cuda_status != cudaSuccess) {
        fprintf(stderr, "rayptr init kernel failed!");
        exit(1);
    }

    int bytesize = sizeof(Float3) * simulation->blocks_per_dim * 2;
    int threads_per_gpublock = 32;
    int gpublocks_needed = ceil((int)((float)NUM_RAYS / (float)threads_per_gpublock));
    int gpublocks_per_SM = ceil((float)gpublocks_needed / (float)n_streams);
    int stride = gpublocks_per_SM * threads_per_gpublock;
    

    printf("Rayinit requires %d bytes of shared mem\n", bytesize);
    printf("Calling initkernel with %d blocks of %d threads per SM. Total: (%d/%d)\n", gpublocks_per_SM, threads_per_gpublock, n_streams*gpublocks_per_SM*threads_per_gpublock, NUM_RAYS);

    Box* boxx;
    boxx = new Box();
    Box* box = simulation->box;
    for (int i = 0; i < n_streams; i++) {
        int offset = stride * i;
        initRayKernel << < gpublocks_per_SM, threads_per_gpublock, bytesize, stream[i] >> > (rayptr, box, focalpoint, offset);
    }
    cudaDeviceSynchronize();

    for (int i = 0; i < n_streams; i++)
        cudaStreamDestroy(stream[i]);

    cuda_status = cudaGetLastError();
    if (cuda_status != cudaSuccess) {
        fprintf(stderr, "rayptr init kernel failed!");
        exit(1);
    }

    

    // FINAL OUTPUT
        
    if (verbose) {
        printf("Ray %d: %f %f %f\n", INDEXA, rayptr[INDEXA].unit_vector.x, rayptr[INDEXA].unit_vector.y, rayptr[INDEXA].unit_vector.z);
        printf("block_indexes: ");
        for (int i = 0; i < MAX_RAY_BLOCKS; i++) {
            printf("%d ", rayptr[INDEXA].block_indexes[i]);
        }
    }

    int duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t0).count();
    printf("\nRaytracer initiated in %d ms\n\n", duration);
}
    
__device__ void colorRay(Ray* ray, uint8_t* image, int index) {

}

__global__ void renderKernel(Ray* rayptr, uint8_t* image, Box* box) {
    int  index = blockIdx.x * blockDim.x + threadIdx.x;
    Ray ray = rayptr[index];
    




    for (int i = 0; i < MAX_RAY_BLOCKS; i++) {
        if (ray.block_indexes[i] == -1)
            break;

        Block* block = &box->blocks[ray.block_indexes[i]];
        if (ray.searchBlock(block, image))
            break;
    }
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