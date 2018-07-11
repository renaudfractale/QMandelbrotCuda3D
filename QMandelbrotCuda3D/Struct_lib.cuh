
#include <iostream> // prompt Output
#include <fstream> //File Output
#include <math.h> //lib mayh
#include <stdio.h> // lib stantard
//#include <cuda_fp16.h> // lib CUDA
#include <windows.h>
#include <fstream>
#include <string>


#define Dim_isFix 0;
#define Dim_end 20.0f;
#define Dim_start -20.0f;
#define Dim_NbStep 4;
#define Dim_step 10.0f;

#define ITER_MAX 255;
#define ITER_MIN 1;
#define ITER_isFix 0;

#define DEV 1;
#define FILTER 0;
#define POWER 2.0f;
#define ISSHOW 1;
#define RMAX 4.0f;

#define NBPOINTS 64;

#define MODEAUTOMANU false;
//Strcture state

 struct struct_Stat_float_T {
	float Xmin = Dim_end
	float Xmax = Dim_start
	float Wmin = Dim_end
	float Wmax = Dim_start
	float Ymin = Dim_end
	float Ymax = Dim_start
	float Zmin = Dim_end
	float Zmax = Dim_start
	float Wstep = 0.0f;
	float Xstep = 0.0f;
	float Ystep = 0.0f;
	float Zstep = 0.0f;
	unsigned long NbPoint = 0;
 };


// struct sur la gestion des dimensions
 struct struct_P_float_T {
	 int isFix = Dim_isFix
	 float start = Dim_start
	 float end = Dim_end
	 int NbStep = Dim_NbStep
	 float step = Dim_step
		int coef = 1;
} ;


struct 	 struct_Iter_T  {
	int max = ITER_MAX
	int min = ITER_MIN
	int isFix = ITER_isFix
} ;

typedef struct 	 {
	char root[100];
	char txt[110];
	char csv[110];
	char stat[110];
	char histo[110];
	char stl[110];
} struct_FileName_T;

// struct sur la gestion paramètres d'entré
struct struct_P_All_T {
	struct_P_float_T X;
	struct_P_float_T Y;
	struct_P_float_T Z;
	struct_P_float_T W;
	struct_Iter_T Iter;
	int dev = DEV
		int filter = FILTER
		float power = POWER
		int isShow = ISSHOW
		struct_FileName_T nameFile;
	float rMax = RMAX
		int NbPointByStep = NBPOINTS
		bool modeAM = MODEAUTOMANU
} ;

struct 	struct_P_Simulation_T {
	//Quaternions
	struct_P_float_T X;
	struct_P_float_T Y;
	struct_P_float_T Z;
	struct_P_float_T W;
	struct_Iter_T Iter;
	float rMax;
	//Parametrer variable systematique
	float power;
	int max;
} ;

struct struct_Q_T {
	float x;
	float y;
	float z;
	float w;
} ;