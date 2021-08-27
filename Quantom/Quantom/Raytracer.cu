#include "Raytracer.cuh"



Ray::Ray(Double3 unit_vector) : unit_vector(unit_vector) {}

__device__ void Ray::findBlockHits(Box* box, Double3 focalpoint) {
    for (int i = 0; i < MAX_RAY_BLOCKS; i++)
		block_indexes[i] = -1;

	int block_index = 0; 

    
	for (int i = 0; i < box->n_blocks; i++) {	// Only works if sorted by y axis first!
		if (hitsBlock(&box->blocks[i], focalpoint)) {
			block_indexes[block_index++] = i;

			if (block_index == MAX_RAY_BLOCKS)
				break;
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

__device__ bool containedBetweenPlanes(double plane_min, double plane_max, double dim_uv) {
    float near = DBL_MIN;
    float far = DBL_MAX;

    float t1 = plane_min / dim_uv;
    float t2 = plane_max / dim_uv;
    float tMin = cudaMin(t1, t2);
    float tMax = cudaMax(t1, t2);
    if (tMin > near) near = tMin;
    if (tMax < far) far = tMax;
    if (near > far || far < 0)
    {
        return false;
    }
    return true;
}

__device__ bool Ray::hitsBlock(Block* block, Double3 focalpoint) {
    
    Double3 blocksize = Double3(BLOCK_LEN, BLOCK_LEN, BLOCK_LEN);
    Double3 min = (block->center - (blocksize * (1/ 2))) - focalpoint;
    Double3 max = (block->center + (blocksize * (1/ 2))) - focalpoint;
    
    float near = DBL_MIN;
    float far = DBL_MAX;

    if (!containedBetweenPlanes(min.x, max.x, unit_vector.x))
        return false;
    if (!containedBetweenPlanes(min.y, max.y, unit_vector.y))
        return false;
    if (!containedBetweenPlanes(min.z, max.z, unit_vector.z))
        return false;


    return true;
}

__global__ void initRayKernel(Ray* rayptr, Box* box, Double3 focalpoint) {
	int  index = blockIdx.x * blockDim.x + threadIdx.x;
    Ray ray = rayptr[index];

	ray.findBlockHits(box, focalpoint);

    rayptr[index] = ray;
}










Raytracer::Raytracer(Simulation* simulation) {
    cuda_status = cudaSetDevice(0);
    if (cuda_status != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!");
        exit(1);
    }



	double box_size = simulation->box_size;
	double principal_point_increment = box_size / RAYS_PER_DIM;

	Ray* host_rayptr = new Ray[NUM_RAYS];
	focalpoint = Double3(0, -box_size / 2 * FOCAL_LEN_RATIO - box_size/2, 0);
	int index = 0;
	for (double z = -box_size / 2 + principal_point_increment / 2; z < box_size / 2; z += principal_point_increment) {
		for (double x = -box_size / 2 + principal_point_increment / 2; x < box_size / 2; x += principal_point_increment) {
			Double3 vector = Double3(x, 0, z) - focalpoint;
			Double3 uv = vector * (1 / vector.len());

            host_rayptr[index++] = Ray(uv);
		}
	}

    cuda_status = cudaMallocManaged(&rayptr, NUM_RAYS * sizeof(Ray));
    cuda_status = cudaMemcpy(rayptr, host_rayptr, NUM_RAYS * sizeof(Ray), cudaMemcpyHostToDevice);


	printf("Allocating %d MB of ram for Rays... ", NUM_RAYS * sizeof(Ray) / 1000000);

    initRayKernel <<< RAYS_PER_DIM, RAYS_PER_DIM, 0>>> (rayptr, simulation->box, focalpoint);
    cudaDeviceSynchronize();
    cuda_status = cudaGetLastError();
    if (cuda_status != cudaSuccess) {
        fprintf(stderr, "rayptr init kernel failed!");
        exit(1);
    }

    printf("Ray 500: %f %f %f\n", rayptr[500].unit_vector.x, rayptr[500].unit_vector.y, rayptr[500].unit_vector.z);
    printf("Rays initiated\n");
}

__global__ void renderKernel(Ray* rayptr, uint8_t* image) {
    int  index = blockIdx.x * 1000 + threadIdx.x;
    Ray ray = rayptr[index];
    Int3 xy_index(index % RAYS_PER_DIM, index / RAYS_PER_DIM, 0);

    double gradient = 255 * ((double)xy_index.x / 1000.);
    
    image[index * 4 + 0] = gradient;
    image[index * 4 + 1] = 0;
    image[index * 4 + 2] = 0;
    image[index * 4 + 3] = 255;
}


uint8_t* Raytracer::render(Simulation* simulation) {
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

    renderKernel << < RAYS_PER_DIM, RAYS_PER_DIM, 0 >> > ( rayptr, cuda_image);
    cudaMemcpy(image, cuda_image, im_bytesize, cudaMemcpyDeviceToHost);


    cudaFree(cuda_image);

    return image;
}