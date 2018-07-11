#include "cuda_runtime.h" //lib W10
#include "device_launch_parameters.h"//lib W10
#include "Struct_lib.cuh"


__device__  void CreateQ_By_float(struct_Q_T *out, float x, float y, float z, float w);
__device__  float  Get_QNorm(struct_Q_T *Q);
__device__ void Get_QPow(struct_Q_T *Q, float pow);