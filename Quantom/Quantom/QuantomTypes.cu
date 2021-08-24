#include "QuantomTypes.cuh"


__host__ __device__ inline double Double3::length() { return (double) sqrtf(x * x + y * y + z * z); }
