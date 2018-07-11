#include "compute.cuh"

// CUDA kernel to Compute itermax of quaternion
__global__ void kernel(const struct_P_Simulation_T *P_Simulation, int *Tab_Iter)
{
	//int Tempindex = 0;
	struct_Q_T Q_Current;
	float w, x, y, z;
	int iter = 0;
	//X
	x = ((float)blockIdx.x)*P_Simulation->X.step + P_Simulation->X.start;
	//Y
	y = ((float)blockIdx.y)*P_Simulation->Y.step + P_Simulation->Y.start;
	//Z
	z = ((float)blockIdx.z)*P_Simulation->Z.step + P_Simulation->Z.start;
	//W
	w = ((float)threadIdx.x)*P_Simulation->W.step + P_Simulation->W.start;

	CreateQ_By_float(&Q_Current, x, y, z, w);

	for (iter = 0; iter <= P_Simulation->Iter.max; iter++)
	{
		Get_QPow(&Q_Current, P_Simulation->power);
		Q_Current.x += x;
		Q_Current.y += y;
		Q_Current.z += z;
		Q_Current.w += w;

		if (Get_QNorm(&Q_Current) > P_Simulation->rMax)
			goto Fin;
	}
Fin:
	if (iter > 0)
		iter--;
	int index = blockIdx.x*P_Simulation->X.coef + blockIdx.y*P_Simulation->Y.coef + blockIdx.z*P_Simulation->Z.coef + threadIdx.x*P_Simulation->W.coef;
	if (index < P_Simulation->max)
		Tab_Iter[index] = iter;// index % 255;
	else
		printf("%d > %d", index, P_Simulation->max);
}



struct_Stat_float_T compute(struct_P_All_T Config, bool state)
{
std::ofstream filetxt;
std::ofstream file;
	
int Tab_Histo[300];
int  Nbpoint_iter = 0;
int NoConfig = 0;
int NbConfig = Config.W.NbStep*Config.X.NbStep*Config.Y.NbStep*Config.Z.NbStep;
//Stat
struct_Stat_float_T Stat;
for (int NoW = 0; NoW < Config.W.NbStep; NoW++)
{
	for (int NoX = 0; NoX < Config.X.NbStep; NoX++)
	{
		for (int NoY = 0; NoY < Config.Y.NbStep; NoY++)
		{
			for (int NoZ = 0; NoZ < Config.Z.NbStep; NoZ++)
			{
				NoConfig++;
				std::cout << "---------------------------------------------------\n";
				std::cout << "Config  " << NoConfig << " sur " << NbConfig << "\n";

				float W = NoW*Config.W.step + Config.W.start;
				float X = NoX*Config.X.step + Config.X.start;
				float Y = NoY*Config.Y.step + Config.Y.start;
				float Z = NoZ*Config.Z.step + Config.Z.start;

				int PasW = Config.NbPointByStep;
				if (Config.W.NbStep == 1)
					PasW = 1;

				int PasX = Config.NbPointByStep;
				if (Config.X.NbStep == 1)
					PasX = 1;

				int PasY = Config.NbPointByStep;
				if (Config.Y.NbStep == 1)
					PasY = 1;

				int PasZ = Config.NbPointByStep;
				if (Config.Z.NbStep == 1)
					PasZ = 1;

				//Taille de tableau
				int max = PasZ*PasY*PasX*PasW;

				if (Config.isShow)
					std::cout << "cudaMallocManaged Config  -->  Start" << "\n";

				struct_P_Simulation_T *P_Simulation;
				int *Tab_Iter;
				// Allocate Unified Memory -- accessible from CPU or GPU
				cudaMallocManaged(&P_Simulation, sizeof(struct_P_Simulation_T));
				cudaMallocManaged(&Tab_Iter, max * sizeof(int));
				if (Config.isShow)
					std::cout << "cudaMallocManaged Config  -->  End " << "\n";

				if (Config.isShow)
					std::cout << "P_Simulation Config  -->  Start" << "\n";
				// Pramatrage de W
				P_Simulation->W.start = W;
				P_Simulation->W.end = W + Config.W.step;
				P_Simulation->W.NbStep = PasW;
				if (PasW == 1)
					P_Simulation->W.step = 0.0f;
				else
					P_Simulation->W.step = (Config.W.step) / (PasW - 1);
				P_Simulation->W.coef = PasX*PasY*PasZ;

				// Pramatrage de X
				P_Simulation->X.start = X;
				P_Simulation->X.end = X + Config.X.step;
				P_Simulation->X.NbStep = PasX;
				if (PasX == 1)
					P_Simulation->X.step = 0.0f;
				else
					P_Simulation->X.step = (Config.X.step) / (PasX - 1);
				P_Simulation->X.coef = PasY*PasZ;

				// Pramatrage de Y
				P_Simulation->Y.start = Y;
				P_Simulation->Y.end = Y + Config.Y.step;
				P_Simulation->Y.NbStep = PasY;
				if (PasY == 1)
					P_Simulation->Y.step = 0.0f;
				else
					P_Simulation->Y.step = (Config.Y.step) / (PasY - 1);
				P_Simulation->Y.coef = PasZ;

				// Pramatrage de Z
				P_Simulation->Z.start = Z;
				P_Simulation->Z.end = Z + Config.Z.step;
				P_Simulation->Z.NbStep = PasZ;
				if (PasZ == 1)
					P_Simulation->Z.step = 0.0f;
				else
					P_Simulation->Z.step = (Config.Z.step) / (PasZ - 1);
				P_Simulation->Z.coef = 1;

				//Stat Step

				Stat.Wstep = P_Simulation->W.step;
				Stat.Xstep = P_Simulation->X.step;
				Stat.Ystep = P_Simulation->Y.step;
				Stat.Zstep = P_Simulation->Z.step;


				//Parametrage Iter
				P_Simulation->Iter.max = Config.Iter.max;
				P_Simulation->Iter.min = Config.Iter.min;

				//Parametrage Power
				P_Simulation->power = Config.power;

				//Parametrage Rmax
				P_Simulation->rMax = Config.rMax;

				//Parametrage max
				P_Simulation->max = max;

				if (Config.isShow)
					std::cout << "P_Simulation Config  -->  End" << "\n";

				if (Config.isShow)
					std::cout << "Tab_Iter and Tab_Histo Init  -->  Start" << "\n";
				for (int i = 0; i < max; i++)
					Tab_Iter[i] = 0;

				for (int i = 0; i <= Config.Iter.max; i++)
					Tab_Histo[i] = 0;
				if (Config.isShow)
					std::cout << "Tab_Iter and Tab_Histo Init -->  End" << "\n";

				if (Config.isShow)
					std::cout << "Compude GPU -->  Start" << "\n";
				dim3 grid(PasX, PasY, PasZ);
				dim3 block(PasW, 1, 1);
				kernel << <grid, block >> >(P_Simulation, Tab_Iter);
				if (Config.isShow)
					std::cout << "Compude GPU -->  End" << "\n";


				if (Config.isShow)
					std::cout << "cudaDeviceSynchronize-->  Start" << "\n";
				cudaDeviceSynchronize();
				if (Config.isShow)
					std::cout << "cudaDeviceSynchronize -->  End" << "\n";

				if (Config.isShow)
					std::cout << "Analyzer Simulation -->  Start" << "\n";
				Nbpoint_iter = 0;
				for (int i = 0; i < max; i++)
				{
					if (Tab_Iter[i] > 0)
						Nbpoint_iter++;
					Tab_Histo[Tab_Iter[i]]++;
				}
				if (Config.isShow)
				{
					std::cout << "Nb point Nbpoint_iter = " << Nbpoint_iter << "\n";
					std::cout << "Soit  :  " << (float)(Nbpoint_iter / ((float)max / 10000.0f)) / 100.0f << "%  soit " << Nbpoint_iter << "pt sur " << max << "pt \n";
					std::cout << "Analyzer Simulation -->  End" << "\n";
				}

				if (Config.isShow)
					std::cout << "Write Histogram -->  Start" << "\n";
				file.open(Config.nameFile.histo, std::ofstream::out | std::ofstream::app);
				file << NoConfig << ";";
				for (int i = 0; i <= Config.Iter.max; i++)
					file << Tab_Histo[i] << ";";
				file << "\n";
				file.close();
				if (Config.isShow)
					std::cout << "Write Histogram -->  End" << "\n";


				file.open(Config.nameFile.csv, std::ofstream::out | std::ofstream::app);
				filetxt.open(Config.nameFile.txt, std::ofstream::out | std::ofstream::app);


				for (int i = 0; i < max; i++)
				{
					int j = i;

					//W
					int iW = 0;
					if (PasW > 1)
					{
						iW = j / P_Simulation->W.coef;
					}
					//printf("index = %d  - Z Tempindex = %d \n", i, Tempindex);
					float w = (float)iW*P_Simulation->W.step + P_Simulation->W.start;
					// on retranche 
					j -= iW*P_Simulation->W.coef;

					//X
					int iX = 0;
					if (PasX > 1)
					{
						iX = j / P_Simulation->X.coef;
					}
					float x = (float)iX*P_Simulation->X.step + P_Simulation->X.start;
					// on retranche 
					j -= iX*P_Simulation->X.coef;

					//Y
					int iY = 0;
					if (PasY > 1)
					{
						iY = j / P_Simulation->Y.coef;
					}
					float y = (float)iY*P_Simulation->Y.step + P_Simulation->Y.start;
					// on retranche 
					j -= iY*P_Simulation->Y.coef;

					//Z
					int iZ = 0;
					if (PasZ > 1)
					{
						iZ = j / P_Simulation->Z.coef;
					}
					float z = (float)iZ*P_Simulation->Z.step + P_Simulation->Z.start;
					// on retranche 
					j -= iZ*P_Simulation->Z.coef;


					int iter = Tab_Iter[i];
					if ((iter >= Config.Iter.min && Config.Iter.isFix == 2) || (iter == Config.Iter.min && Config.Iter.isFix == 1))
					{
						int filter = Config.filter;
						if (state == true)
						{
							if (FilterQ_H(&filter, &iX, &iY, &iZ, &iW, iter, P_Simulation))
							{
								file << x << ";" << y << ";" << z << ";" << w << ";" << iter << ";\n";
								int NbDim = 0;
								if (PasW > 1)
								{
									filetxt << w << ";";
									NbDim++;
								}

								if (PasX > 1)
								{
									filetxt << x << ";";
									NbDim++;
								}

								if (PasY > 1)
								{
									filetxt << y << ";";
									NbDim++;
								}

								if (PasZ > 1)
								{
									filetxt << z << ";";
									NbDim++;
								}

								if (NbDim >= 3)
									filetxt << iter << "\n";
								else
									filetxt << ((float)iter) / ((float)Config.Iter.max)*3.0f << "\n";
							}
						}
						if (w > Stat.Wmax)
							Stat.Wmax = w;
						if (w < Stat.Wmin)
							Stat.Wmin = w;

						if (x > Stat.Xmax)
							Stat.Xmax = x;
						if (x < Stat.Xmin)
							Stat.Xmin = x;

						if (y > Stat.Ymax)
							Stat.Ymax = y;
						if (y < Stat.Ymin)
							Stat.Ymin = y;

						if (z > Stat.Zmax)
							Stat.Zmax = z;
						if (z < Stat.Zmin)
							Stat.Zmin = z;

						Stat.NbPoint++;
					}

				}
				file.close();
				filetxt.close();

				if (Config.isShow)
					std::cout << "Clear Mem + Reste  -->  Start" << "\n";
				cudaFree(P_Simulation);
				cudaFree(Tab_Iter);
				//cudaDeviceReset();
				if (Config.isShow)
					std::cout << "Clear Mem + Reste  -->  End" << "\n";
			}
		}
	}
}

return Stat;

}