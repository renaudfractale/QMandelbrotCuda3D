#include "cuda_runtime.h" //lib W10
#include "device_launch_parameters.h"//lib W10
#include "device_lib.cuh"



//__managed__  struct_P_Simulation_T *P_Simulation;
//__managed__  int *Tab_Iter;

__host__  void CreateQ_By_float_H(struct_Q_T *out, float x, float y, float z, float w);
__host__  float  Get_QNorm_H(struct_Q_T *Q);
__host__ void Get_QPow_H(struct_Q_T *Q, float pow);
__host__ int  GetQIter_H(struct_P_Simulation_T *P_Simulation_DEVICE, int  *x_filter, int  *y_filter, int *z_filter, int *w_filter);
__host__ bool  FilterQ_H(int *Filter, int *Nx, int *Ny, int *Nz, int *Nw, int iter, struct_P_Simulation_T *P_Simulation);