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
		double d1 = (box->blocks[block_indexes[i - 1]].center - focalpoint).length();
		double d2 = (box->blocks[block_indexes[i]].center - focalpoint).length();

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

__device__ bool Ray::hitsBlock(Block* block, Double3 focalpoint) {
    Double3 blocksize = Double3(block->length, block->length, block->length);
    Double3 min = (block->center - (blocksize * (1/ 2))) - focalpoint;
    Double3 max = (block->center + (blocksize * (1/ 2))) - focalpoint;
    float near = DBL_MIN;
    float far = DBL_MAX;

    // X
    float t1 = min.x / unit_vector.x;
    float t2 = max.x / unit_vector.x;
    float tMin = cudaMin(t1, t2);
    float tMax = cudaMax(t1, t2);
    if (tMin > near) near = tMin;
    if (tMax < far) far = tMax;
    if (near > far || far < 0)
    {
        return false;
    }

    // Y
    t1 = min.y / unit_vector.y;
    t2 = max.y / unit_vector.y;
    tMin = cudaMin(t1, t2);
    tMax = cudaMax(t1, t2);
    if (tMin > near) near = tMin;
    if (tMax < far) far = tMax;
    if (near > far || far < 0)
    {
        return false;
    }

    // Z
    t1 = min.z / unit_vector.z;
    t2 = max.z / unit_vector.z;
    tMin = cudaMin(t1, t2);
    tMax = cudaMax(t1, t2);
    if (tMin > near) near = tMin;
    if (tMax < far) far = tMax;
    if (near > far || far < 0)
    {
        return false;
    }

    //point1 = origin + direction * near;
    //point2 = origin + direction * far;
    return true;
}

__global__ void initRayKernel(Ray* rayptr, Box* box, Double3 focalpoint) {
	int  index = blockIdx.x * THREADS_PER_BLOCK + threadIdx.x;
	Ray ray = rayptr[index];

	ray.findBlockHits(box, focalpoint);
	
}










Raytracer::Raytracer(Simulation* simulation) {
	double box_size = simulation->box_size;

	double principal_point_increment = box_size / RAYS_PER_DIM;

	rayptr = new Ray[RAYS_PER_DIM * RAYS_PER_DIM];

	Double3 focal_point(0, box_size / 2 * 4, 0);
	int index = 0;
	for (int z = -box_size / 2 + principal_point_increment / 2; z < box_size / 2; z += principal_point_increment) {
		for (int x = -box_size / 2 + principal_point_increment / 2; x < box_size / 2; x += principal_point_increment) {
			Double3 vector = Double3(x, 0, z) - focal_point;
			Double3 uv = vector * (1 / vector.length());

			rayptr[index++] = Ray(uv);
		}
	}



	//cudaMallocManaged(&rayptr, NUM_RAYS * sizeof(Ray));
	printf("Allocating %d MB of ram for Rays... ", NUM_RAYS * sizeof(Ray) / 1000000);

    //initRayKernel <<< 1, THREADS_PER_BLOCK, 0>>> (rayptr, simulation->box, focal_point);

    printf("Rays initiated\n");
}

