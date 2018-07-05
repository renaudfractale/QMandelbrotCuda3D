#include "cuda_runtime.h" //lib W10
#include "device_launch_parameters.h"//lib W10
#include <iostream> // prompt Output
#include <fstream> //File Output
#include <math.h> //lib mayh
#include <stdio.h> // lib stantard
#include <cuda_fp16.h> // lib CUDA


// Pour X,Y,Z et W
typedef struct 	struct_P_float {
	float start;
	float end;
	int NbPoints;
	float step;
} struct_P_float_T;


typedef struct 	struct_P_Power {
	float value;
} struct_P_Power_T;

typedef struct 	struct_P_Iter {
	int start;
	int end;
} struct_P_Iter_T;


typedef struct 	struct_P_Rlimit {
	float value;
} struct_P_Rlimit_T;

typedef struct 	struct_P_Simulation {
	//Quaternions
	struct_P_Power_T X;
	struct_P_float_T Y;
	struct_P_float_T Z;
	struct_P_float_T W;
	//Parametre Fixe
	struct_P_Iter_T Iter;
	float Rlimit;
	//Parametrer variable systematique
	float Power;
} struct_P_Simulation_T;

typedef struct 	struct_Q {
	float x;
	float y;
	float z;
	float w;
} struct_Q_T;

struct_P_Simulation_T *P_Simulation_DEVICE;
short *Tab_Iter;


__device__  void CreateQ_By_float(struct_Q_T *out, float x, float y, float z, float w)
{
	out->x = x;
	out->y = y;
	out->z = z;
	out->w = w;
}

__device__  float  Get_QNorm(struct_Q_T *Q)
{
	return sqrtf(Q->x*Q->x + Q->y*Q->y + Q->z*Q->z + Q->w*Q->w);
}

__device__ void Get_QPow(struct_Q_T *Q, float pow)
{
	float A = Get_QNorm(Q);
	float theta = 0.0f;
	float B = 0.0f;
	float R = 0.0f;
	if (pow > 0.0f && A>0.000001f)
	{
		float coef = 1.0f;
		if (A<1.0f)
		{
			//printf("%f *******\n", A);
			coef = 1 / A;
			Q->x *= coef;
			Q->y *= coef;
			Q->z *= coef;
			Q->z *= coef;

		}
		A = Get_QNorm(Q);
		//printf("%f +++++++++\n", A);
		theta = acosf(Q->w / A)*pow;
		B = sqrt(A*A - Q->w*Q->w);
		R = exp2f(logf(A / coef)* pow);
		Q->x = R*sinf(theta)*(Q->x / B);
		Q->y = R*sinf(theta)*(Q->y / B);
		Q->z = R*sinf(theta)*(Q->z / B);
		Q->z = R*cosf(theta);

	}
	else
	{
		//printf("%f --------\n", A);
		Q->w = 0.0f;
		Q->x = 0.0f;
		Q->y = 0.0f;
		Q->z = 0.0f;

	}
}

//
//__device__  void  Get_QAdd(struct_Q_T *Q1, struct_Q_T *Q2)
//{
//	Q1->x += Q2->x;
//	Q1->y += Q2->y;
//	Q1->z += Q2->z;
//	Q1->w += Q2->w;
//}





// CUDA kernel to Compute itermax of quaternion
__global__ void kernel(const struct_P_Simulation_T *P_Simulation, short *Tab_Iter)
{
	//int Tempindex = 0;
	struct_Q_T Q_Current;
	float w, x, y, z;
	int iter = 0;
	//X
	x = P_Simulation->X.value;
	//Y
	y = ((float)blockIdx.y)*P_Simulation->Y.step + P_Simulation->Y.start;
	//Z
	z = ((float)blockIdx.x)*P_Simulation->Z.step + P_Simulation->Z.start;
	//W
	w = ((float)threadIdx.x)*P_Simulation->W.step + P_Simulation->W.start;

	CreateQ_By_float(&Q_Current, x, y, z, w);

	for (iter = 0; iter <= P_Simulation->Iter.end; iter++)
	{
		Get_QPow(&Q_Current, P_Simulation->Power);
		Q_Current.x += x;
		Q_Current.y += y;
		Q_Current.z += z;
		Q_Current.w += w;

		if (Get_QNorm(&Q_Current) > P_Simulation->Rlimit)
			goto Fin;
	}
Fin:
	if (iter > 0)
		iter--;
	Tab_Iter[(blockIdx.y * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x] = (short)iter;
}

int main(int argc, char *argv[])
{

	//#################### CONFIG par DEFAUT #################
	struct_P_float_T ParameterDelaults;
	// ******************** NbPoints **************
	ParameterDelaults.NbPoints = 10;
	char Str_NbPoints[] = "-NbPoints";
	//********************* Parameter y,z,w *************
	ParameterDelaults.start = -3.0f;
	char Str_Q_start[] = "-Q_start";
	ParameterDelaults.end = 3.0f;
	char Str_Q_end[] = "-Q_end";
	ParameterDelaults.step = (ParameterDelaults.end - ParameterDelaults.start) / ((float)ParameterDelaults.NbPoints - 1);
	//char Str_Q_step[] = "-Q_step";
	//*********************  NameFile Output *************
	char NameFile[110];
	char NameFile_csv[110];
	char NameFile_histo[110];
	strcpy(NameFile, "OutputFile");
	char Str_NameFile[] = "-o";
	//*********************  X value *************
	float X = -0.3375f;
	char Str_X[] = "-X";

	float POWER = 2.0;
	char Str_Power[] = "-Power";

	char Str_H[] = "-h";
	char Str_Help[] = "--help";

	int itermax = 255;
	float Rmax = 4.0;

	int dev = 0;
	char Str_dev[] = "-dev";

	bool IsShow = 0;
	char Str_IsShow[] = "-IsShow";

	if (argc > 1)
	{
		if (argc % 2 == 0)
		{
			for (int i = 1; i <= argc; i++)
			{
				if (strcmp(argv[i], Str_Help) == 0 || strcmp(argv[i], Str_H) == 0)
				{
					int count;
					cudaGetDeviceCount(&count);

					std::cout << "Help :  \n";
					std::cout << "      -NbPoints : numbers of points \n";
					std::cout << "                  Ctrt 01 : if start = end , NbPoints must be equal to 1 \n";
					std::cout << "                  Ctrt 02 : if start < end , NbPoints must be sup to 1 \n";
					std::cout << "                  Ctrt 03 : NbPoints must be type int\n";
					std::cout << "      -Q_start : value  of start \n";
					std::cout << "                  Ctrt 03 :  start >= end \n";
					std::cout << "                  Ctrt 03 : start must be type float\n";
					std::cout << "      -Q_end : value  of end \n";
					std::cout << "                  Ctrt 03 :  start >= end \n";
					std::cout << "                  Ctrt 03 : end must be type float\n";
					std::cout << "      -o : Output File \n";
					std::cout << "                  Ctrt 05 :  len  must be inf 100 char \n";
					std::cout << "      -X : value of x \n";
					std::cout << "                  Ctrt 03 : x must be type float\n";
					std::cout << "      -Power : value of Power \n";
					std::cout << "                  Ctrt 03 : Power must be type float\n";
					std::cout << "                  Ctrt 03 : Power must be sup 0.0 \n";
					std::cout << "      -dev : index of device GPU, 0 by default \n";
					std::cout << "            --> you have " << count << " device(s) in your PC \n";
					std::cout << "                  Ctrt 03 : dev must be type int\n";
					std::cout << "                  Ctrt 03 : dev must be sup or equal 0 \n";
					std::cout << "      -IsShow : bool if show message in cmd, true (1) by default \n";
					std::cout << "                  Ctrt 03 : IsShow must be type int (0 or 1)\n";
					std::cout << "                  Ctrt 03 : IsShow must be sup or equal 0 \n";
					std::cout << "      -h / --help : show help \n";
					std::cout << "      Example :\n";
					std::cout << "               Programme.exe -X 0.3375 -Q_start -3.0 -Q_end 3.0 -NbPoints 4 -Power 2.5 -o FileOutput \n";
					std::cout << "               Programme.exe -X 0.3375 -Q_start -4.0 -Q_end 4.0 -NbPoints 4 -o FileOutput2 -IsShow 0 -dev 1 -Power 2.0\n";
					std::cout << "      Version : 0.1 du 04 Juillet 2018\n";
					std::cout << "      Auteur : Renaud HENRY\n";
					std::cout << "      siteweb : http://fractale.io/ \n";
					return 0;
				}
			}
			std::cout << "Error 00 : Argument impaire" << "\n";
			return -1;
		}
		for (int i = 1; i < argc; i += 2)
		{
			std::cout << "Analyse du couple d'arguments :  " << argv[i] << " " << argv[i + 1] << "\n";
			if (strcmp(argv[i], Str_Q_start) == 0)
			{
				ParameterDelaults.start = (float)atof(argv[i + 1]);
				if (errno)
				{
					std::cout << "Error 02 " << Str_Q_start << ": value is not float " << "\n";
					return -1;
				}
			}
			else if ((strcmp(argv[i], Str_Q_end) == 0))
			{
				ParameterDelaults.end = (float)atof(argv[i + 1]);
				if (errno)
				{
					std::cout << "Error 02 " << Str_Q_end << ": value is not type float " << "\n";
					return -1;
				}
			}
			else if ((strcmp(argv[i], Str_X) == 0))
			{
				X = (float)atof(argv[i + 1]);
				if (errno)
				{
					std::cout << "Error 02 " << Str_X << ": value is not type float " << "\n";
					return -1;
				}
			}
			else if ((strcmp(argv[i], Str_NbPoints) == 0))
			{
				ParameterDelaults.NbPoints = atoi(argv[i + 1]);
				if (errno)
				{
					std::cout << "Error 02 " << Str_NbPoints << ": value is not type int " << "\n";
					return -1;
				}
			}
			else if ((strcmp(argv[i], Str_Power) == 0))
			{
				POWER = (float)atof(argv[i + 1]);
				if (errno)
				{
					std::cout << "Error 09 " << Str_Power << ": value is not type float " << "\n";
					return -1;
				}
			}
			else if ((strcmp(argv[i], Str_IsShow) == 0))
			{
				int temp = atoi(argv[i + 1]);
				if (errno)
				{
					std::cout << "Error 12 " << Str_IsShow << ": value is not type int " << "\n";
					return -1;
				}
				if (temp == 0)
					IsShow = false;
				else if (temp == 1)
					IsShow = true;
				else
				{
					std::cout << "Error 13 " << Str_IsShow << ": must be equal at 1 ou 0" << "\n";
					return -1;
				}
			}
			else if ((strcmp(argv[i], Str_dev) == 0))
			{
				dev = atoi(argv[i + 1]);
				if (errno)
				{
					std::cout << "Error 11 " << Str_dev << ": value is not type int " << "\n";
					return -1;
				}
				int count;
				cudaGetDeviceCount(&count);
				if (count < dev)
				{
					std::cout << "Error 12 " << Str_dev << ": value of dev > countdevice, " << dev << " > " << count << "\n";
					return -1;
				}

			}
			else if ((strcmp(argv[i], Str_NameFile) == 0))
			{
				if (strlen(argv[i + 1])<100)
					strcpy(NameFile, argv[i + 1]);
				else
				{
					std::cout << "Error 07 strlen fileOutput must be inf to 100 signe \n";
					return -1;
				}
			}
			else
			{
				std::cout << "Warning 08  Arg not know : " << argv[i] << " " << argv[i + 1] << "\n";
			}
		}
	}
	if (POWER <= 0.0f)
	{
		std::cout << "Error 10 POWER < 0.0f : " << POWER << "<" << 0.0f << "\n";
		return -1;
	}
	if (ParameterDelaults.end < ParameterDelaults.start)
	{
		std::cout << "Error 03 end < start : " << ParameterDelaults.end << "<" << ParameterDelaults.start << "\n";
		return -1;
	}
	if (ParameterDelaults.end == ParameterDelaults.start && ParameterDelaults.NbPoints != 1)
	{
		std::cout << "Warning  04 end == start  and NbPoints !=1:  So NbPoints force to 1\n";
		ParameterDelaults.NbPoints = 1;
	}
	if (ParameterDelaults.NbPoints > 1)
	{
		ParameterDelaults.step = (ParameterDelaults.end - ParameterDelaults.start) / ((float)ParameterDelaults.NbPoints);
		if (ParameterDelaults.step < 0.0001f)
		{
			std::cout << "Error 05 step < 0.0001 : " << ParameterDelaults.end << "<" << 0.0001f << "\n";
			std::cout << "step must be sup at 0.0001 \n";
			return -1;
		}
	}
	else if (ParameterDelaults.end != ParameterDelaults.start)
	{
		std::cout << "Error 06 NbPoints < 1:  So NbPoints must be sup 0 \n";
		return -1;
	}
	/****** Print Parmaters Used ******/
	std::cout << "Parameters Current : " << "\n";
	std::cout << "		Q_start = " << ParameterDelaults.start << ", Q_end = " << ParameterDelaults.end << ", Q_Step = " << ParameterDelaults.step << ", Nbpoints = " << ParameterDelaults.NbPoints << "\n";
	std::cout << "		x_value = " << X << "\n";
	std::cout << "		FileOutput = " << NameFile << "\n";
	std::cout << "		IterMax_Max = " << itermax << ", Rmax = " << Rmax << "\n";

	/********  Clear File ************/
	std::ofstream file;
	strcpy(NameFile_csv, NameFile);
	strcat(NameFile_csv, ".csv");
	file.open(NameFile_csv);
	file << "X;Y;Z;W;iter;\n";
	file.close();

	strcpy(NameFile_histo, NameFile);
	strcat(NameFile_histo, ".histo");
	file.open(NameFile_histo);
	file << "index;";
	for (int i = 0; i <= itermax; i++)
		file << i << ";";
	file << "\n";
	file.close();

	//#################### Constante(s) #################

	const int  maxMaster = ParameterDelaults.NbPoints * ParameterDelaults.NbPoints  * ParameterDelaults.NbPoints;
	const int NbPoints = 128;
	const int  maxMinor = NbPoints * NbPoints  * NbPoints;

	//#################### Variables(s) #################

	int Y = 0;
	int Z = 0;
	int W = 0;

	int Tab_Histo[300];
	int  Nbpoint_iter = 0;

	cudaSetDevice(dev);
	for (int index = 0; index < maxMaster; index++)
	{
		int indexTemp = index;

		Y = indexTemp / (ParameterDelaults.NbPoints * ParameterDelaults.NbPoints);
		indexTemp = indexTemp - (ParameterDelaults.NbPoints * ParameterDelaults.NbPoints)*Y;

		Z = indexTemp / (ParameterDelaults.NbPoints);
		indexTemp = indexTemp - (ParameterDelaults.NbPoints * Z);

		W = indexTemp;


		if (IsShow)
			std::cout << "cudaMallocManaged Config  -->  Start" << "\n";
		// Allocate Unified Memory -- accessible from CPU or GPU
		cudaMallocManaged(&P_Simulation_DEVICE, sizeof(struct_P_Simulation_T));
		cudaMallocManaged(&Tab_Iter, maxMinor * sizeof(short));
		if (IsShow)
			std::cout << "cudaMallocManaged Config  -->  End " << "\n";
		if (IsShow)
			std::cout << "P_Simulation Config  -->  Start" << "\n";
		//Parametrage Iter
		P_Simulation_DEVICE->Iter.end = itermax;
		P_Simulation_DEVICE->Iter.start = 10;

		//Parametrage Power
		P_Simulation_DEVICE->Power = POWER;

		//Parametrage Rmax
		P_Simulation_DEVICE->Rlimit = Rmax;

		//Parametrage X
		P_Simulation_DEVICE->X.value = X;

		//Parametrage Y
		P_Simulation_DEVICE->Y.start = (float)Y*ParameterDelaults.step + ParameterDelaults.start;
		P_Simulation_DEVICE->Y.end = (float)(Y + 1)*ParameterDelaults.step + ParameterDelaults.start;
		P_Simulation_DEVICE->Y.NbPoints = NbPoints;
		P_Simulation_DEVICE->Y.step = (P_Simulation_DEVICE->Y.end - P_Simulation_DEVICE->Y.start) / (P_Simulation_DEVICE->Y.NbPoints - 1);

		//Parametrage Z
		P_Simulation_DEVICE->Z.start = (float)Z*ParameterDelaults.step + ParameterDelaults.start;
		P_Simulation_DEVICE->Z.end = (float)(Z + 1)*ParameterDelaults.step + ParameterDelaults.start;
		P_Simulation_DEVICE->Z.NbPoints = NbPoints;
		P_Simulation_DEVICE->Z.step = (P_Simulation_DEVICE->Z.end - P_Simulation_DEVICE->Z.start) / (P_Simulation_DEVICE->Z.NbPoints - 1);

		//Parametrage W
		P_Simulation_DEVICE->W.start = (float)W*ParameterDelaults.step + ParameterDelaults.start;
		P_Simulation_DEVICE->W.end = (float)(W + 1)*ParameterDelaults.step + ParameterDelaults.start;
		P_Simulation_DEVICE->W.NbPoints = NbPoints;
		P_Simulation_DEVICE->W.step = (P_Simulation_DEVICE->W.end - P_Simulation_DEVICE->W.start) / (P_Simulation_DEVICE->W.NbPoints - 1);
		if (IsShow)
			std::cout << "P_Simulation Config  -->  End" << "\n";

		std::cout << "P_Simulation Config No " << index + 1 << " sur  " << maxMaster << "\n";

		if (IsShow)
			std::cout << "Tab_Iter and Tab_Histo Init  -->  Start" << "\n";
		for (int i = 0; i < maxMinor; i++) {
			Tab_Iter[i] = (short)0;

		}

		for (int i = 0; i <= itermax; i++)
			Tab_Histo[i] = 0;
		if (IsShow)
			std::cout << "Tab_Iter and Tab_Histo Init -->  End" << "\n";

		if (IsShow)
			std::cout << "Compude GPU -->  Start" << "\n";
		int NbThreadPerBlock = NbPoints;
		int NbBlockPerGrid = NbPoints;
		dim3 grid(NbBlockPerGrid, NbBlockPerGrid, 1);
		dim3 block(NbThreadPerBlock, 1, 1);
		kernel << <grid, block >> >(P_Simulation_DEVICE, Tab_Iter);
		if (IsShow)
			std::cout << "Compude GPU -->  End" << "\n";

		if (IsShow)
			std::cout << "cudaDeviceSynchronize-->  Start" << "\n";
		cudaDeviceSynchronize();
		if (IsShow)
			std::cout << "cudaDeviceSynchronize -->  End" << "\n";

		if (IsShow)
			std::cout << "Analyzer Simulation -->  Start" << "\n";
		Nbpoint_iter = 0;
		for (int i = 0; i < maxMinor; i++)
		{
			if (Tab_Iter[i] > 0)
				Nbpoint_iter++;
			Tab_Histo[Tab_Iter[i]]++;
		}
		if (IsShow)
			std::cout << "Nb point Nbpoint_iter = " << Nbpoint_iter << "\n";
		if (IsShow)
			std::cout << "Soit  :  " << (float)(Nbpoint_iter / (maxMinor / 10000)) / 100.0f << "%  soit " << Nbpoint_iter << "pt sur " << maxMinor << "pt \n";
		if (IsShow)
			std::cout << "Analyzer Simulation -->  End" << "\n";

		if (IsShow)
			std::cout << "Write Histogram -->  Start" << "\n";
		file.open(NameFile_histo, std::ofstream::out | std::ofstream::app);
		file << index << ";";
		for (int i = 0; i <= itermax; i++)
			file << Tab_Histo[i] << ";";
		file << "\n";
		file.close();
		if (IsShow)
			std::cout << "Write Histogram -->  End" << "\n";

		if (IsShow)
			std::cout << "Write csv -->  Start" << "\n";
		file.open(NameFile_csv, std::ofstream::out | std::ofstream::app);

		for (int i = 0; i < maxMinor; i++)
		{
			int j = i;
			//X
			float x = P_Simulation_DEVICE->X.value;

			//Y
			float y = ((float)(j / (NbPoints*NbPoints))*P_Simulation_DEVICE->Y.step) + P_Simulation_DEVICE->Y.start;
			// on retranche 
			j -= (j / (NbPoints*NbPoints))*(NbPoints*NbPoints);

			//Z
			float z = (float)(j / NbPoints)*P_Simulation_DEVICE->Z.step + P_Simulation_DEVICE->Z.start;
			// on retranche Q2
			j -= (j / NbPoints)*NbPoints;

			//W
			//printf("index = %d  - Z Tempindex = %d \n", i, Tempindex);
			float w = (((float)j)*P_Simulation_DEVICE->W.step) + P_Simulation_DEVICE->W.start;

			short iter = Tab_Iter[i];
			if (iter > 0)
			{
				file << x << ";" << y << ";" << z << ";" << w << ";" << iter << ";\n";
			}

		}
		file.close();
		if (IsShow)
			std::cout << "Write csv -->  End" << "\n";

		if (IsShow)
			std::cout << "Clear Mem + Reste  -->  Start" << "\n";
		cudaFree(P_Simulation_DEVICE);
		cudaFree(Tab_Iter);
		//cudaDeviceReset();
		if (IsShow)
			std::cout << "Clear Mem + Reste  -->  End" << "\n";
	}
	return 0;
}