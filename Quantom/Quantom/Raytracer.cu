#include "Raytracer.cuh"



Ray::Ray(Float3 unit_vector, Float3 origin) : unit_vector(unit_vector), origin(origin) {}

__device__ void Ray::findBlockHits(Box* box, Float3 focalpoint) {
    for (int i = 0; i < MAX_RAY_BLOCKS; i++)
		block_indexes[i] = -1;

    int block_index = 0;


	// Only works if sorted by y axis first!                      

    

    

    for (int y = 0; y < box->blocks_per_dim; y++) {
        for (int z = 0; z < box->blocks_per_dim; z++) {
            for (int x = 0; x < box->blocks_per_dim; x++) {
                if (block_index == MAX_RAY_BLOCKS)
                    break;

                int index = z * box->blocks_per_dim * box->blocks_per_dim + y * box->blocks_per_dim + x;


                if (hitsBlock(&box->blocks[index], focalpoint)) {
                    block_indexes[block_index] = index;

                    block_index++;
                    
                        
                }

            }
        }
    }


    
	for (int i = 1; i < block_index; i++) {		// I think only 1 swap will be necessary...
		float d1 = (box->blocks[block_indexes[i - 1]].center - focalpoint).len();
		float d2 = (box->blocks[block_indexes[i]].center - focalpoint).len();

		if (d2 < d1) {
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

    //Float3 min = block->center - Float3(BLOCK_LEN / 2, BLOCK_LEN / 2, BLOCK_LEN / 2);
    //Float3 max = block->center + Float3(BLOCK_LEN / 2, BLOCK_LEN / 2, BLOCK_LEN / 2);

    Float3 offset(BODY_RADIUS, BODY_RADIUS, BODY_RADIUS);   // Radius of largest molecule!
    //Float3 offset(0, 0, 0);
    Float3 min = block->center - Float3(FOCUS_LEN_HALF, FOCUS_LEN_HALF, FOCUS_LEN_HALF) - offset *2;    //I have no idea why we need *2
    Float3 max = block->center + Float3(FOCUS_LEN_HALF, FOCUS_LEN_HALF, FOCUS_LEN_HALF) + offset *2;

    float tmin = -DBL_MAX;
    float tmax = DBL_MAX;

    /*if (blockIdx.x * 1000 + threadIdx.x == INDEXA) {
        printf("boxmin: %.2f %.2f %.2f\n", min.x, min.y, min.z);
        printf("boxmax: %.2f %.2f %.2f\n", max.x, max.y, max.z);
        //printf("Tmin: %f \t tmax: %f \n", tmin, tmax);
    }*/

    for (int dim = 0; dim < 3; dim++) {
        float invD = unit_vector.at(dim) != 0 ? 1.0f / unit_vector.at(dim) : 999999999;
        //float invD = 1.0f / unit_vector.at(dim);
        float t0 = (min.at(dim) - focalpoint.at(dim)) * invD;
        float t1 = (max.at(dim) - focalpoint.at(dim)) * invD;
        //if (blockIdx.x * 1000 + threadIdx.x == INDEXA)
          //  printf("t0 = %.4f\t t1 = %.4f\n", t0, t1);
        if (invD < 0.0f) {
            float temp = t1;
            t1 = t0;
            t0 = temp;
        }

        tmin = t0 > tmin ? t0 : tmin;
        tmax = t1 < tmax ? t1 : tmax;
        //if (blockIdx.x * 1000 + threadIdx.x == INDEXA)
          //  printf("%f %f\n", tmin, tmax);


        if (tmax <= tmin) {    // was <=
            return false;
        }
    } 
    /*if (blockIdx.x * 1000 + threadIdx.x == INDEXA) {
        printf("HIT! Block body count: %d. First body pos : %.6f %.6f %.6f\n\n\n", block->n_bodies, block->bodies[0].pos.x, block->bodies[0].pos.y, block->bodies[0].pos.z);
    }*/
        
    return true;
    /*
    if (!containedBetweenPlanes(min.x, max.x, unit_vector.x, focalpoint.x))
        return false;
    if (!containedBetweenPlanes(min.y, max.y, unit_vector.y, focalpoint.y))
        return false;
    if (!containedBetweenPlanes(min.z, max.z, unit_vector.z, focalpoint.z))
        return false;


    return true;*/
}
    

__device__ bool Ray::hitsBody(SimBody* body) {
    if (distToPoint(body->pos) < BODY_RADIUS)
        return true;
    return false;
}

__device__ bool Ray::moleculeCollisionHandling(SimBody* body, MoleculeLibrary* mol_library, uint8_t* image) {
    Molecule* mol = &mol_library->molecules[0];

    const int infinity = 9999999;

    float closest_collision = infinity;
    Atom closest_atom;
    closest_atom.pos = Float3(0, infinity, 0);  // Make sure its infinitely far away in y direction.




    Float3 molecule_tilt_vector = Float3(0, 1, 0).rotateAroundOrigin(body->rotation);
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

    if (closest_atom.pos.y != infinity) {
        int index = blockIdx.x * blockDim.x + threadIdx.x;

        image[index * 4 + 0] = closest_atom.color[0];
        image[index * 4 + 1] = closest_atom.color[1];
        image[index * 4 + 2] = closest_atom.color[2];
        image[index * 4 + 3] = 255;
        return true;
    }
    
    
    return false;
}

__device__ float Ray::distToSphereIntersect(Atom* atom) {
    Float3 projection_on_ray = origin + unit_vector * ((atom->pos - origin).dot(unit_vector) / unit_vector.dot(unit_vector));
    float center_to_projection = (projection_on_ray - atom->pos).len();
    float projection_to_intersect = sqrtf(atom->radius * atom->radius - center_to_projection * center_to_projection);
    return (projection_on_ray - origin).len() - projection_to_intersect;
}

__global__ void initRayKernel(Ray* rayptr, Box* box, Float3 focalpoint) {
	int  index = blockIdx.x * blockDim.x + threadIdx.x;
    Ray ray = rayptr[index];

	ray.findBlockHits(box, focalpoint);

    rayptr[index] = ray;
}


    

        
        
        








    





Raytracer::Raytracer(Simulation* simulation) {
    setGPU();



	float box_size = simulation->box_size;
    float base = -(box_size + FOCUS_LEN) / 2.f;
	float principal_point_increment = (box_size + FOCUS_LEN) / (float)RAYS_PER_DIM;

	Ray* host_rayptr = new Ray[NUM_RAYS];
	focalpoint = Float3(0, -(box_size  / 2.f) * FOCAL_LEN_RATIO - box_size, 0);  

    printf("\n\nFocal point: %.3f %.3f %.3f\n", focalpoint.x, focalpoint.y, focalpoint.z);
	int index = 0;
    for (int z_index = 0; z_index < RAYS_PER_DIM; z_index++) {
        for (int x_index = 0; x_index < RAYS_PER_DIM; x_index++) {

            float z = base + principal_point_increment * (float)z_index;
            float x = base + principal_point_increment * (float)x_index;
			Float3 vector = Float3(x, 0, z) - focalpoint;
			Float3 uv = vector * (1.f / vector.len());     
            host_rayptr[index++] = Ray(uv, focalpoint);
		}
	}

    cuda_status = cudaMallocManaged(&rayptr, NUM_RAYS * sizeof(Ray));
    cuda_status = cudaMemcpy(rayptr, host_rayptr, NUM_RAYS * sizeof(Ray), cudaMemcpyHostToDevice);


	printf("Allocating %d MB of ram for Rays... \n", NUM_RAYS * sizeof(Ray) / 1000000);

    initRayKernel <<< RAYS_PER_DIM, RAYS_PER_DIM>>> (rayptr, simulation->box, focalpoint);
    cudaDeviceSynchronize();

    cuda_status = cudaGetLastError();
    if (cuda_status != cudaSuccess) {
        fprintf(stderr, "rayptr init kernel failed!");
        exit(1);
    }
        
        
    int indexa = INDEXA;
    printf("Ray %d: %f %f %f\n", indexa, rayptr[indexa].unit_vector.x, rayptr[indexa].unit_vector.y, rayptr[indexa].unit_vector.z);
    printf("block_indexes: ");
    for (int i = 0; i < MAX_RAY_BLOCKS; i++) {
        printf("%d ", rayptr[indexa].block_indexes[i]);
    }


    printf("Rays initiated\n\n");
}
    
__device__ void colorRay(Ray* ray, uint8_t* image, int index) {

}

__global__ void renderKernel(Ray* rayptr, uint8_t* image, Box* box, MoleculeLibrary* mol_library) {
    int  index = blockIdx.x * blockDim.x + threadIdx.x;
    Ray ray = rayptr[index];
    


    for (int i = 0; i < MAX_RAY_BLOCKS; i++) {
        if (ray.block_indexes[i] == -1)
            break;

        Block* block = &box->blocks[ray.block_indexes[i]];
        for (int j = 0; j < MAX_FOCUS_BODIES; j++) {

            if (block->focus_bodies[j].molecule_type == UNUSED_BODY) {  // We can do this because bodies are loaded from index 0 at each timestep. MIGHT NEED TO CHANGE IN THE FUTURE!
                break;
            }

            if (ray.hitsBody(&block->focus_bodies[j])) {
                if (ray.moleculeCollisionHandling(&block->focus_bodies[j], mol_library, image)) {
                    return;
                }
            }
        }
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


    uint8_t* image = new uint8_t[NUM_RAYS * 4];
    uint8_t* cuda_image;
    int im_bytesize = NUM_RAYS * 4 * sizeof(uint8_t);
    cudaMallocManaged(&cuda_image, im_bytesize);
    //printf("Allocating %d KB of ram for image... ", im_bytesize / 1000);
    Float3 a(1, 0, 1);
    renderKernel << < RAYS_PER_DIM, RAYS_PER_DIM, 0>>> ( rayptr, cuda_image, simulation->box, simulation->mol_library);
    cudaMemcpy(image, cuda_image, im_bytesize, cudaMemcpyDeviceToHost);
    //exit(1);

    cudaFree(cuda_image);
    cudaStreamDestroy(renderstream);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    printf("\rRender time: %4d ms  ", duration.count());
    // First render: 666 ms

    return image;
}               