#include "Raytracer.cuh"



Ray::Ray(Double3 unit_vector, Double3 origin) : unit_vector(unit_vector), origin(origin) {}

__device__ void Ray::findBlockHits(Box* box, Double3 focalpoint) {
    for (int i = 0; i < MAX_RAY_BLOCKS; i++)
		block_indexes[i] = -1;

    int block_index = 0;


	// Only works if sorted by y axis first!                      

    

    int bpd = BOX_LEN_CUDA / BLOCK_LEN_CUDA;

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
		double d1 = (box->blocks[block_indexes[i - 1]].center - focalpoint).len();
		double d2 = (box->blocks[block_indexes[i]].center - focalpoint).len();

		if (d2 < d1) {
			int tmp = block_indexes[i - 1];
			block_indexes[i - 1] = block_indexes[i];
			block_indexes[i] = tmp;
		}
	}
      
}

__device__ double cudaMax(double a, double b) { 
    if (a > b)
        return a;
    return b;
}
__device__ double cudaMin(double a, double b) {
    if (a < b)
        return a;
    return b;
}

__device__ bool containedBetweenPlanes(double plane_min, double plane_max, double dir_dim, double origin_dim) {
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

__device__ bool Ray::hitsBlock(Block* block, Double3 focalpoint) {
    double a = BLOCK_LEN_CUDA;               // FUCK YOU VISUAL STUDIO
    Double3 blocksize = Double3(a, a, a);

    Double3 min = block->center - Double3(BLOCK_LEN_CUDA / 2, BLOCK_LEN_CUDA / 2, BLOCK_LEN_CUDA / 2);
    Double3 max = block->center + Double3(BLOCK_LEN_CUDA / 2, BLOCK_LEN_CUDA / 2, BLOCK_LEN_CUDA / 2);

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
    if (distToPoint(body->pos) < 0.005)
        return true;
    return false;
}

__global__ void initRayKernel(Ray* rayptr, Box* box, Double3 focalpoint) {
	int  index = blockIdx.x * blockDim.x + threadIdx.x;
    Ray ray = rayptr[index];

	ray.findBlockHits(box, focalpoint);

    rayptr[index] = ray;
}





















Raytracer::Raytracer(Simulation* simulation) {
    setGPU();



	double box_size = simulation->box_size;
	double principal_point_increment = box_size / (double)RAYS_PER_DIM;

	Ray* host_rayptr = new Ray[NUM_RAYS];
	focalpoint = Double3(0, -(box_size / 2.f) * FOCAL_LEN_RATIO - box_size, 0);  

    printf("\n\nFocal point: %.3f %.3f %.3f\n", focalpoint.x, focalpoint.y, focalpoint.z);
	int index = 0;
    for (int z_index = 0; z_index < RAYS_PER_DIM; z_index++) {
        for (int x_index = 0; x_index < RAYS_PER_DIM; x_index++) {

            double z = -box_size / 2.f + principal_point_increment * (double)z_index;
            double x = -box_size / 2.f + principal_point_increment * (double)x_index;
			Double3 vector = Double3(x, 0, z) - focalpoint;

           

			Double3 uv = vector * (1.f / vector.len());
            if (index == INDEXA) {
                printf("%.2f %.2f %.2f\n", vector.x, vector.y, vector.z);
                printf("%.2f %.2f %.2f\n", uv.x, uv.y, uv.z);
            }
                
            host_rayptr[index++] = Ray(uv, focalpoint);
		}
	}

    cuda_status = cudaMallocManaged(&rayptr, NUM_RAYS * sizeof(Ray));
    cuda_status = cudaMemcpy(rayptr, host_rayptr, NUM_RAYS * sizeof(Ray), cudaMemcpyHostToDevice);


	printf("Allocating %d MB of ram for Rays... \n", NUM_RAYS * sizeof(Ray) / 1000000);

    initRayKernel <<< RAYS_PER_DIM, RAYS_PER_DIM, 0>>> (rayptr, simulation->box, focalpoint);

    cuda_status = cudaGetLastError();
    if (cuda_status != cudaSuccess) {
        fprintf(stderr, "rayptr init kernel failed!");
        exit(1);
    }
    cudaDeviceSynchronize();
       

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

__global__ void renderKernel(Ray* rayptr, uint8_t* image, Box* box, Double3 focalpoint) {
    int  index = blockIdx.x * blockDim.x + threadIdx.x;

    Ray ray = rayptr[index];
    //Int3 xy_index(blockIdx.x, threadIdx.x, 0);
     


    for (int i = 0; i < MAX_RAY_BLOCKS; i++) {
        if (ray.block_indexes[i] == -1)
            break;

        Block* block = &box->blocks[ray.block_indexes[i]];

        for (int j = 0; j < block->n_bodies; j++) {
            
            if (ray.hitsBody(&block->bodies[j])) {

                image[index * 4 + 0] = 220;
                image[index * 4 + 1] = 0;
                image[index * 4 + 2] = 0;
                image[index * 4 + 3] = 255;
                return;
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

    uint8_t* image = new uint8_t[NUM_RAYS * 4];
    uint8_t* cuda_image;
    int im_bytesize = NUM_RAYS * 4 * sizeof(uint8_t);
    cudaMallocManaged(&cuda_image, im_bytesize);
    //printf("Allocating %d KB of ram for image... ", im_bytesize / 1000);
    Double3 a(1, 0, 1);
    renderKernel << < RAYS_PER_DIM, RAYS_PER_DIM, 0 >> > ( rayptr, cuda_image, simulation->box, focalpoint);
    cudaMemcpy(image, cuda_image, im_bytesize, cudaMemcpyDeviceToHost);


    cudaFree(cuda_image);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    printf("Render time: %d\r", duration.count());
    // First render: 666 ms

    return image;
}