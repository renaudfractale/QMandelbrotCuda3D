#include "cuda_runtime.h" //lib W10
#include "device_launch_parameters.h"//lib W10
#include "compute.cuh"


int main(int argc, char *argv[])
{
	struct_P_All_T Config;
    int state= menu( argc, argv, &Config);

	if (state != 0)
		return state;

	/********  Clear File ************/
	std::ofstream filetxt;
	filetxt.open(Config.nameFile.txt);
	filetxt.close();

	std::ofstream file;
	file.open(Config.nameFile.csv);
	file << "X;Y;Z;W;iter;\n";
	file.close();


	file.open(Config.nameFile.histo);
	file << "index;";
	for (int i = 0; i <= Config.Iter.max; i++)
		file << i << ";";
	file << "\n";
	file.close();

	//Affiche Stat
	FILE *fileman;
	char line[1000];
	fileman = fopen(Config.nameFile.stat, "r");
	if (fileman != NULL)
	{
		while (std::fgets(line, 1000, fileman))
		{
			std::cout << line;
		}
		fclose(fileman);
	}
	else
		std::cout << "Error  : "<< Config.nameFile.stat <<" no found" << "\n";
	cudaSetDevice(Config.dev);


	struct_Stat_float_T Stat;

	if(Config.modeAM==false)
		Stat =compute(Config, true);
	else
	{
		bool Dim[] = {false,false ,false ,false };
		
		if (Config.W.isFix == 2)
			Dim[0] = true;
		if (Config.X.isFix == 2)
			Dim[1] = true;
		if (Config.Y.isFix == 2)
			Dim[2] = true;
		if (Config.Z.isFix == 2)
			Dim[3] = true;

		int NbDim = 0;

		for (int i = 0; i <= 3; i++)
		{
			if (Dim[i]==true)
			{
				NbDim++;
				switch (i)
				{
				case 0:
					Config.W.isFix = 2;
					Config.W.start = -20.0f;
					Config.W.end = 20.0f;
					Config.W.NbStep = 2;
					Config.W.step = (Config.W.end - Config.W.start) / ((float)Config.W.NbStep);
					break;
				case 1:
					Config.X.isFix = 2;
					Config.X.start = -20.0f;
					Config.X.end = 20.0f;
					Config.X.NbStep = 2;
					Config.X.step = (Config.X.end - Config.X.start) / ((float)Config.X.NbStep);
					break;
				case 2:
					Config.Y.isFix = 2;
					Config.Y.start = -20.0f;
					Config.Y.end = 20.0f;
					Config.Y.NbStep = 2;
					Config.Y.step = (Config.Y.end - Config.Y.start) / ((float)Config.Y.NbStep);
					break;
				case 3:
					Config.Z.isFix = 2;
					Config.Z.start = -20.0f;
					Config.Z.end = 20.0f;
					Config.Z.NbStep = 2;
					Config.Z.step = (Config.Z.end - Config.Z.start) / ((float)Config.Z.NbStep);
					break;
				default:
					break;
				}
			}	
		}
		//si et selement si il y a 3 dimentions
		if (NbDim == 3)
		{
			float pas = 10;
			float OldPas = 20;
			int nb_iter = 10;
			for (int iter = 0; iter <= nb_iter; iter++)
			{
				Stat = compute(Config,false);
				// rétro action
				for (int i = 0; i <= 3; i++)
				{
					if (Dim[i] == true)
					{
						switch (i)
						{
						case 0:
							if (Stat.Wmin == Config.W.start)
							{
								Config.W.start -= OldPas;
								if (nb_iter != i)
									Config.W.start += pas;
							}
							else if(nb_iter!=i)
								Config.W.start +=pas;

							if (Stat.Wmax == Config.W.end)
							{
								Config.W.end += OldPas;
								if (nb_iter != i)
									Config.W.end -= pas;
							}								
							else if(nb_iter != i)
								Config.W.end -= pas;
							Config.W.step = (Config.W.end - Config.W.start) / ((float)Config.W.NbStep);
							break;
						case 1:
							if (Stat.Xmin == Config.X.start)
							{
								Config.X.start -= OldPas;
								if (nb_iter != i)
									Config.X.start += pas;
							}
							else if (nb_iter != i)
								Config.X.start += pas;
							if (Stat.Xmax == Config.X.end)
							{
								Config.X.end += OldPas;
								if (nb_iter != i)
									Config.X.end -= pas;
							}
							else if (nb_iter != i)
								Config.X.end -= pas;
							Config.X.step = (Config.X.end - Config.X.start) / ((float)Config.X.NbStep);
							break;
						case 2:
							if (Stat.Ymin == Config.Y.start)
							{
								Config.Y.start -= OldPas;
								if (nb_iter != i)
									Config.Y.start += pas;
							}
							else if (nb_iter != i)
								Config.Y.start += pas;
							if (Stat.Ymax == Config.Y.end)
							{
								Config.Y.end += OldPas;
								if (nb_iter != i)
									Config.Y.end -= pas;
							}
							else if (nb_iter != i)
								Config.Y.end -= pas;
							Config.Y.step = (Config.Y.end - Config.Y.start) / ((float)Config.Y.NbStep);
							break;
						case 3:
							if (Stat.Zmin == Config.Z.start)
							{
								Config.Z.start -= OldPas;
								if (nb_iter != i)
									Config.Z.start += pas;
							}
							else if (nb_iter != i)
								Config.Z.start += pas;
							if (Stat.Zmax == Config.Z.end)
							{
								Config.Z.end += OldPas;
								if (nb_iter != i)
									Config.Z.end -= pas;
							}
							else if (nb_iter != i)
								Config.Z.end -= pas;
							Config.Z.step = (Config.Z.end - Config.Z.start) / ((float)Config.Z.NbStep);
							break;
						default:
							break;
						}
					}

				}

				file.open(Config.nameFile.stat, std::ofstream::out | std::ofstream::app);
				file << "Statistiques : \n";
				file << "				X min = " << Stat.Xmin << "\n";
				file << "				X max = " << Stat.Xmax << "\n";
				file << "				Y min = " << Stat.Ymin << "\n";
				file << "				Y max = " << Stat.Ymax << "\n";
				file << "				Z min = " << Stat.Zmin << "\n";
				file << "				Z max = " << Stat.Zmax << "\n";
				file << "				W min = " << Stat.Wmin << "\n";
				file << "				W max = " << Stat.Wmax << "\n";
				file << "				X step = " << Stat.Xstep << "\n";
				file << "				Y step = " << Stat.Ystep << "\n";
				file << "				Z step = " << Stat.Zstep << "\n";
				file << "				W step = " << Stat.Wstep << "\n";
				file << "				NbPoint plot = " << Stat.NbPoint << "\n";
				file.close();


				OldPas = pas;
				pas = pas / 2.0f;
			}

			for (int i = 0; i <= 3; i++)
			{
				if (Dim[i] == true)
				{
					switch (i)
					{
					case 0:
						Config.W.NbStep = 5;
						Config.W.step = (Config.W.end - Config.W.start) / ((float)Config.W.NbStep);
						break;
					case 1:
						Config.X.NbStep = 5;
						Config.X.step = (Config.X.end - Config.X.start) / ((float)Config.X.NbStep);
						break;
					case 2:
						Config.Y.NbStep = 5;
						Config.Y.step = (Config.Y.end - Config.Y.start) / ((float)Config.Y.NbStep);
						break;
					case 3:
						Config.Z.NbStep = 5;
						Config.Z.step = (Config.Z.end - Config.Z.start) / ((float)Config.Z.NbStep);
						break;
					default:
						break;
					}
				}

			}
			if(Stat.NbPoint!=0)
				Stat = compute(Config, true);
		}
		else
		{
			std::cout << "Error  mode: mode auto only 3 dim\n";
			return -1;
		}

	}
	
	file.open(Config.nameFile.stat, std::ofstream::out | std::ofstream::app);
	file << "Statistiques : \n";
	file << "				X min = " << Stat.Xmin << "\n";
	file << "				X max = " << Stat.Xmax << "\n";
	file << "				Y min = " << Stat.Ymin << "\n";
	file << "				Y max = " << Stat.Ymax << "\n";
	file << "				Z min = " << Stat.Zmin << "\n";
	file << "				Z max = " << Stat.Zmax << "\n";
	file << "				W min = " << Stat.Wmin << "\n";
	file << "				W max = " << Stat.Wmax << "\n";
	file << "				X step = " << Stat.Xstep << "\n";
	file << "				Y step = " << Stat.Ystep << "\n";
	file << "				Z step = " << Stat.Zstep << "\n";
	file << "				W step = " << Stat.Wstep << "\n";
	file << "				NbPoint plot = " << Stat.NbPoint << "\n";
	file.close();




	return 0;//Fin du programme
}