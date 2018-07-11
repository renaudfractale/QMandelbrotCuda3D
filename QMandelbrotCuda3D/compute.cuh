#include "cuda_runtime.h" //lib W10
#include "device_launch_parameters.h"//lib W10
#include "menu_lib.cuh"


__global__ void kernel(const struct_P_Simulation_T *P_Simulation, int *Tab_Iter);
struct_Stat_float_T compute(struct_P_All_T Config, bool state);
