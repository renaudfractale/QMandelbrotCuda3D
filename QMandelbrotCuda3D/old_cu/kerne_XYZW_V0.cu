#include "cuda_runtime.h" //lib W10
#include "device_launch_parameters.h"//lib W10
#include <iostream> // prompt Output
#include <fstream> //File Output
#include <math.h> //lib mayh
#include <stdio.h> // lib stantard
#include <cuda_fp16.h> // lib CUDA
#include <windows.h>
#include <fstream>
#include <string>
#define Dim_isFix 0;
#define Dim_end 20.0f;
#define Dim_start -20.0f;
#define Dim_NbStep 4;
#define Dim_step 10.0f;

#define ITER_MAX 255;
#define ITER_MIN 10;
#define ITER_isFix 0;

#define DEV 1;
#define FILTER 0;
#define POWER 2.0f;
#define ISSHOW 1;
#define RMAX 4.0f;

#define NBPOINTS 64;

// struct sur la gestion des dimensions
typedef struct 	struct_P_float {
	int isFix = Dim_isFix
	float start = Dim_start
	float end = Dim_end
	int NbStep = Dim_NbStep
	float step = Dim_step
} struct_P_float_T;

typedef struct 	struct_Iter {
	int max = ITER_MAX
	int min = ITER_MIN
	int isFix = ITER_isFix
} struct_Iter_T;

typedef struct 	struct_FileName {
	char root[100];
	char txt[110];
	char csv[110];
	char stat[110];
	char histo[110];
	char stl[110];
} struct_FileName_T;

// struct sur la gestion paramètres d'entré
typedef struct 	struct_P_All {
	struct_P_float_T X;
	struct_P_float_T Y;
	struct_P_float_T Z;
	struct_P_float_T W;
	struct_Iter_T Iter;
	int dev = DEV
	int filter=FILTER
	float power = POWER
	int isShow = ISSHOW
	struct_FileName_T nameFile;
	float rMax = RMAX
	int NbPointByStep = NBPOINTS
} struct_P_All_T;


int main(int argc, char *argv[])
{
	//Config
		struct_P_All_T Config;
		strcpy(Config.nameFile.root, "O");
	//Arg Help
		char Str_H[] = "-h";
		char Str_Help[] = "--help";
	
	//Arg X
		char Str_xFix[] = "-x";
		char Str_xMax[] = "-xmax";
		char Str_xMin[] = "-xmin";
		char Str_xNbStep[] = "-xNbStep";
	//Arg W
		char Str_wFix[] = "-w";
		char Str_wMax[] = "-wmax";
		char Str_wMin[] = "-wmin";
		char Str_wNbStep[] = "-wNbStep";
	//Arg Y
		char Str_yFix[] = "-y";
		char Str_yMax[] = "-ymax";
		char Str_yMin[] = "-ymin";
		char Str_yNbStep[] = "-yNbStep";
	//Arg Z
		char Str_zFix[] = "-z";
		char Str_zMax[] = "-zmax";
		char Str_zMin[] = "-zmin";
		char Str_zNbStep[] = "-zNbStep";
	//Arg Iter
		char Str_IterFix[] = "-iter";
		char Str_IterMax[] = "-iterMax";
		char Str_IterMin[] = "-iterMin";
	//Arg Dev
		char Str_dev[] = "-device";
	//Arg Filter
		char Str_filter[] = "-filter";
	//Arg Power
		char Str_power[] = "-power";
	//Arg IsShow
		char Str_IsShow[] = "-isShow";
	//Arg Output File
		char Str_Out[] = "-o";
	//Arg Rmax
		char Str_rMax[] = "-rMax";
	//Si il y a des Arguments
	if (argc > 1)
	{
		// dédection la commande -help pu -h
		for (int i = 1; i <= argc; i++)
		{
			if (strcmp(argv[i], Str_Help) == 0 || strcmp(argv[i], Str_H) == 0)
			{
				//Affiche Help
				FILE *fileman;
				char line[1000];
				fileman = fopen("man", "r");
				if (fileman != NULL)
				{
					while (std::fgets(line, 1000, fileman))
					{
						std::cout << line;
					}
					fclose(fileman);
					return 0; //Fin du programme
				}
				else
				{
					std::cout << "Error  : man no found" << "\n";
					return -1;
				}			
			}
		}
		std::cout << "¨Pass 1" << "\n";
		//Verification : si nb d'arguments est paire --> erreur
		if (argc % 2 == 0)
		{
			std::cout << "Error 00 : Argument impaire" << "\n";
			return -1;
		}
		else //: si nb d'arguments est impaire --> fonctionement normale
		{
			std::cout << "¨Pass 2" << "\n";
			for (int i = 1; i < argc; i += 2)
			{
				std::cout << "Analyse du couple d'arguments :  " << argv[i] << " " << argv[i + 1] << "\n";
				if (strcmp(argv[i], Str_xFix) == 0)
				{
					if (Config.X.isFix == 0) //Si premier config
					{
						Config.X.isFix = 1;
						float value = (float)atof(argv[i + 1]);
						if (errno)
						{
							std::cout << "Error 01 " << Str_xFix << ": value is not float " << "\n";
							return -1;
						}
						Config.X.start = value;
						Config.X.end = value;
						Config.X.NbStep = 1;
						Config.X.step = 1.0f;
					}
					else
					{
						std::cout << "Error 02 " << Str_xFix << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_wFix) == 0)
				{
					if (Config.W.isFix == 0) //Si premier config
					{
						Config.W.isFix = 1;
						float value = (float)atof(argv[i + 1]);
						if (errno)
						{
							std::cout << "Error 01 " << Str_wFix << ": value is not float " << "\n";
							return -1;
						}
						Config.W.start = value;
						Config.W.end = value;
						Config.W.NbStep = 1;
						Config.W.step = 1.0f;
					}
					else
					{
						std::cout << "Error 02 " << Str_wFix << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_yFix) == 0)
				{
					if (Config.Y.isFix == 0) //Si premier config
					{
						Config.Y.isFix = 1;
						float value = (float)atof(argv[i + 1]);
						if (errno)
						{
							std::cout << "Error 01 " << Str_yFix << ": value is not float " << "\n";
							return -1;
						}
						Config.Y.start = value;
						Config.Y.end = value;
						Config.Y.NbStep = 1;
						Config.Y.step = 1.0f;
					}
					else
					{
						std::cout << "Error 02 " << Str_yFix << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_zFix) == 0)
				{
					if (Config.Z.isFix == 0) //Si premier config
					{
						Config.Z.isFix = 1;
						float value = (float)atof(argv[i + 1]);
						if (errno)
						{
							std::cout << "Error 01 " << Str_zFix << ": value is not float " << "\n";
							return -1;
						}
						Config.Z.start = value;
						Config.Z.end = value;
						Config.Z.NbStep = 1;
						Config.Z.step = 1.0f;
					}
					else
					{
						std::cout << "Error 02 " << Str_zFix << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_xMax) == 0)
				{
					if (Config.X.isFix == 0 || Config.X.isFix == 2)
					{
						Config.X.isFix = 2;
						float value = (float)atof(argv[i + 1]);
						if (errno)
						{
							std::cout << "Error 01 " << Str_xMax << ": value is not float " << "\n";
							return -1;
						}
						Config.X.end = value;
					}
					else
					{
						std::cout << "Error 02 " << Str_xMax << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_wMax) == 0)
				{
					if (Config.W.isFix == 0 || Config.W.isFix == 2)
					{
						Config.W.isFix = 2;
						float value = (float)atof(argv[i + 1]);
						if (errno)
						{
							std::cout << "Error 01 " << Str_wMax << ": value is not float " << "\n";
							return -1;
						}
						Config.W.end = value;
					}
					else
					{
						std::cout << "Error 02 " << Str_wMax << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_yMax) == 0)
				{
					if (Config.Y.isFix == 0 || Config.Y.isFix == 2)
					{
						Config.Y.isFix = 2;
						float value = (float)atof(argv[i + 1]);
						if (errno)
						{
							std::cout << "Error 01 " << Str_yMax << ": value is not float " << "\n";
							return -1;
						}
						Config.Y.end = value;
					}
					else
					{
						std::cout << "Error 02 " << Str_yMax << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_zMax) == 0)
				{
					if (Config.Z.isFix == 0 || Config.Z.isFix == 2)
					{
						Config.Z.isFix = 2;
						float value = (float)atof(argv[i + 1]);
						if (errno)
						{
							std::cout << "Error 01 " << Str_zMax << ": value is not float " << "\n";
							return -1;
						}
						Config.Z.end = value;
					}
					else
					{
						std::cout << "Error 02 " << Str_zMax << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_xMin) == 0)
				{
					if (Config.X.isFix == 0 || Config.X.isFix == 2)
					{
						Config.X.isFix = 2;
						float value = (float)atof(argv[i + 1]);
						if (errno)
						{
							std::cout << "Error 01 " << Str_xMin << ": value is not float " << "\n";
							return -1;
						}
						Config.X.start = value;
					}
					else
					{
						std::cout << "Error 02 " << Str_xMin << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_wMin) == 0)
				{
					if (Config.W.isFix == 0 || Config.W.isFix == 2)
					{
						Config.W.isFix = 2;
						float value = (float)atof(argv[i + 1]);
						if (errno)
						{
							std::cout << "Error 01 " << Str_wMin << ": value is not float " << "\n";
							return -1;
						}
						Config.W.start = value;
					}
					else
					{
						std::cout << "Error 02 " << Str_wMin << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_yMin) == 0)
				{
					if (Config.Y.isFix == 0 || Config.Y.isFix == 2)
					{
						Config.Y.isFix = 2;
						float value = (float)atof(argv[i + 1]);
						if (errno)
						{
							std::cout << "Error 01 " << Str_yMin << ": value is not float " << "\n";
							return -1;
						}
						Config.Y.start = value;
					}
					else
					{
						std::cout << "Error 02 " << Str_yMin << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_zMin) == 0)
				{
					if (Config.Z.isFix == 0 || Config.Z.isFix == 2)
					{
						Config.Z.isFix = 2;
						float value = (float)atof(argv[i + 1]);
						if (errno)
						{
							std::cout << "Error 01 " << Str_zMin << ": value is not float " << "\n";
							return -1;
						}
						Config.Z.start = value;
					}
					else
					{
						std::cout << "Error 02 " << Str_zMin << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_xNbStep) == 0)
				{
					if (Config.X.isFix == 0 || Config.X.isFix == 2)
					{
						Config.X.isFix = 2;
						int value = atoi(argv[i + 1]);
						if (errno || value<=0)
						{
							std::cout << "Error 03 " << Str_xNbStep << ": value is not int or value <= 0 " << "\n";
							return -1;
						}
						Config.X.NbStep = value;
					}
					else
					{
						std::cout << "Error 02 " << Str_xNbStep << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_wNbStep) == 0)
				{
					if (Config.W.isFix == 0 || Config.W.isFix == 2)
					{
						Config.W.isFix = 2;
						int value = atoi(argv[i + 1]);
						if (errno || value<=0)
						{
							std::cout << "Error 03 " << Str_wNbStep << ": value is not int or value <= 0 " << "\n";
							return -1;
						}
						Config.W.NbStep = value;
					}
					else
					{
						std::cout << "Error 02 " << Str_wNbStep << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_yNbStep) == 0)
				{
					if (Config.Y.isFix == 0 || Config.Y.isFix == 2)
					{
						Config.Y.isFix = 2;
						int value = atoi(argv[i + 1]);
						if (errno || value<=0)
						{
							std::cout << "Error 03 " << Str_yNbStep << ": value is not int or value <= 0 " << "\n";
							return -1;
						}
						Config.Y.NbStep = value;
					}
					else
					{
						std::cout << "Error 02 " << Str_yNbStep << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_zNbStep) == 0)
				{
					if (Config.Z.isFix == 0 || Config.Z.isFix == 2)
					{
						Config.Z.isFix = 2;
						int value = atoi(argv[i + 1]);
						if (errno || value<=0)
						{
							std::cout << "Error 03 " << Str_yNbStep << ": value is not int or value <= 0 " << "\n";
							return -1;
						}
						Config.Z.NbStep = value;
					}
					else
					{
						std::cout << "Error 02 " << Str_yNbStep << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_IterFix) == 0)
				{
					if (Config.Iter.isFix == 0) //Si premier config
					{
						Config.Iter.isFix = 1;
						int value = atoi(argv[i + 1]);
						if (errno || value<=0)
						{
							std::cout << "Error 03 " << Str_IterFix << ": value is not int  or value <= 0" << "\n";
							return -1;
						}
						Config.Iter.max = value;
						Config.Iter.min = value;
					}
					else
					{
						std::cout << "Error 02 " << Str_IterFix << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_IterMax) == 0)
				{
					if (Config.Iter.isFix == 0 || Config.Iter.isFix == 2) //Si premier config
					{
						Config.Iter.isFix = 2;
						int value = atoi(argv[i + 1]);
						if (errno || value <= 0)
						{
							std::cout << "Error 03 " << Str_IterMax << ": value is not int  or value <= 0" << "\n";
							return -1;
						}
						Config.Iter.max = value;
					}
					else
					{
						std::cout << "Error 02 " << Str_IterMax << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_IterMin) == 0)
				{
					if (Config.Iter.isFix == 0 || Config.Iter.isFix == 2) //Si premier config
					{
						Config.Iter.isFix = 2;
						int value = atoi(argv[i + 1]);
						if (errno || value <= 0)
						{
							std::cout << "Error 03 " << Str_IterMin << ": value is not int  or value <= 0" << "\n";
							return -1;
						}
						Config.Iter.min = value;
					}
					else
					{
						std::cout << "Error 02 " << Str_IterMin << ": value already assigned " << "\n";
						return -1;
					}
				}
				else if (strcmp(argv[i], Str_dev) == 0)
				{
					int value = atoi(argv[i + 1]);
					if (errno)
					{
						std::cout << "Error 01 " << Str_dev << ": value is not int" << "\n";
						return -1;
					}
					Config.dev = value;
				}
				else if (strcmp(argv[i], Str_filter) == 0)
				{
					int value = atoi(argv[i + 1]);
					if (errno)
					{
						std::cout << "Error 01 " << Str_filter << ": value is not int" << "\n";
						return -1;
					}
					Config.filter = value;
				}
				else if (strcmp(argv[i], Str_power) == 0)
				{
					float value = (float)atof(argv[i + 1]);
					if (errno)
					{
						std::cout << "Error 01 " << Str_filter << ": value is not float" << "\n";
						return -1;
					}
					Config.power = value;
				}
				else if (strcmp(argv[i], Str_IsShow) == 0)
				{
					int value = atoi(argv[i + 1]);
					if (errno)
					{
						std::cout << "Error 01 " << Str_IsShow << ": value is not int" << "\n";
						return -1;
					}
					Config.isShow = value;
				}
				else if (strcmp(argv[i], Str_rMax) == 0)
				{
					float value = (float)atof(argv[i + 1]);
					if (errno)
					{
						std::cout << "Error 01 " << Str_rMax << ": value is not float" << "\n";
						return -1;
					}
					Config.rMax = value;
				}
				else if (strcmp(argv[i], Str_Out) == 0)
				{
					if (strlen(argv[i + 1])<100)
						strcpy(Config.nameFile.root, argv[i + 1]);
					else
					{
						std::cout << "Error 04 strlen fileOutput must be inf to 100 signe \n";
						return -1;
					}
				}
				else
				{
					std::cout << "Warning 05  Arg not know : " << argv[i] << " " << argv[i + 1] << "\n";
				}
			}
		}
	}
	bool IsErrors = false;
	// Verification W
	if (Config.W.isFix == 1 && Config.W.start == Config.W.end)
	{
		;//OK
	}
	else if(Config.W.isFix == 2 || Config.W.isFix == 0)
	{
		Config.W.isFix = 2;
		if (Config.W.start < Config.W.end && Config.W.NbStep>1)
		{
			Config.W.step = (Config.W.end - Config.W.start) / ((float)Config.W.NbStep);
		}
		else
		{
			std::cout << "Error W  :  wmax must be sup wmin AND nbStep must be sup 1\n";
			IsErrors=true;
		}
	}
	else
	{
		std::cout << "Error W  :  Error unknow\n";
		IsErrors=true;
	}
	// Verification X
	if (Config.X.isFix == 1 && Config.X.start == Config.X.end)
	{
		;//OK
	}
	else if (Config.X.isFix == 2 || Config.X.isFix == 0)
	{
		Config.X.isFix = 2; 
		if (Config.X.start < Config.X.end && Config.X.NbStep>1)
		{
			Config.X.step = (Config.X.end - Config.X.start) / ((float)Config.X.NbStep);
		}
		else
		{
			std::cout << "Error X  :  xmax must be sup xmin AND nbStep must be sup 1\n";
			IsErrors=true;
		}
	}
	else
	{
		std::cout << "Error X  :  Error unknow\n";
		IsErrors=true;
	}
	// Verification Y
	if (Config.Y.isFix == 1 && Config.Y.start == Config.Y.end)
	{
		;//OK
	}
	else if (Config.Y.isFix == 2 || Config.Y.isFix == 0)
	{
		Config.Y.isFix = 2; 
		if (Config.Y.start < Config.Y.end && Config.Y.NbStep>1)
		{
			Config.Y.step = (Config.Y.end - Config.Y.start) / ((float)Config.Y.NbStep);
		}
		else
		{
			std::cout << "Error Y  :  ymax must be sup ymin AND nbStep must be sup 1\n";
			IsErrors=true;
		}
	}
	else
	{
		std::cout << "Error Y  :  Error unknow\n";
		IsErrors=true;
	}
	// Verification Z
	if (Config.Z.isFix == 1 && Config.Z.start == Config.Z.end)
	{
		;//OK
	}
	else if (Config.Z.isFix == 2 || Config.Z.isFix == 0)
	{
		Config.Z.isFix = 2; 
		if (Config.Z.start < Config.Z.end && Config.Z.NbStep>1)
		{
			Config.Z.step = (Config.Z.end - Config.Z.start) / ((float)Config.Z.NbStep);
		}
		else
		{
			std::cout << "Error Z  :  zmax must be sup zmin AND nbStep must be sup 1\n";
			IsErrors=true;
		}
	}
	else
	{
		std::cout << "Error Z  :  Error unknow\n";
		IsErrors=true;
	}
	// Verification Iter
	if (Config.Iter.isFix == 1)
	{
		if (Config.Iter.min == Config.Iter.max && Config.Iter.max > 0)
		{
			; //OK
		}
		else
		{
			std::cout << "Error Iter  :  value must be sup at 0\n";
			IsErrors=true;
		}
	}
	else if (Config.Iter.isFix == 2 || Config.Iter.isFix == 0)
	{
		Config.Iter.isFix = 2;
		if (Config.Iter.min < Config.Iter.max)
		{
			; // OK
		}
		else
		{
			std::cout << "Error Iter  :  max must be sup min\n";
			IsErrors=true;
		}
	}
	else
	{
		std::cout << "Error Iter  :  Error unknow\n";
		IsErrors=true;
	}
	// Vérification dev
	int count;
	cudaGetDeviceCount(&count);
	if (Config.dev >= 0 && Config.dev < count)
	{
		; //OK
	}
	else
	{
		std::cout << "Error dev :  dev must be between 0 and "<< count -1 <<"\n";
		IsErrors=true;
	}
	//Verification filter
	if (Config.filter >= 0 && Config.filter <= 2)
	{
		; //OK
	}
	else
	{
		std::cout << "Error filter :  filter must be between 0 and 2\n";
		IsErrors=true;
	}
	//Verification power
	if (Config.power >= 2.0f && Config.power <= 50.0f)
	{
		; //OK
	}
	else
	{
		std::cout << "Error power :  power must be between 2.0 and 50.0\n";
		IsErrors=true;
	}
	//Verification IsShow
	if (Config.isShow == 0 || Config.isShow == 1)
	{
		; //OK
	}
	else
	{
		std::cout << "Error isShow :  isShow must be between 0 and 1 \n";
		IsErrors=true;
	}
	//Verification rMax
	if (Config.rMax >0.0f)
	{
		; //OK
	}
	else
	{
		std::cout << "Error rMax :  rMax must be sup 0.0 \n";
		IsErrors = true;
	}
	if (IsErrors)
		return -1;

	// creation des fichiers :
	strcpy(Config.nameFile.csv, Config.nameFile.root);
	strcpy(Config.nameFile.histo, Config.nameFile.root);
	strcpy(Config.nameFile.stat, Config.nameFile.root);
	strcpy(Config.nameFile.stl, Config.nameFile.root);
	strcpy(Config.nameFile.txt, Config.nameFile.root);
	
	strcat(Config.nameFile.csv, ".csv");
	strcat(Config.nameFile.histo, ".histo");
	strcat(Config.nameFile.stat, ".stat");
	strcat(Config.nameFile.stl, ".stl");
	strcat(Config.nameFile.txt, ".txt");

	//Affichage de la config
	std::ofstream file;
	file.open(Config.nameFile.stat);
	file << "Parameters Current : " << "\n";
	file << "				W_start = " << Config.W.start << ", W_end = " << Config.W.end << ", W_Step = " << Config.W.step << ", W_NbStep = " << Config.W.NbStep << "\n";
	file << "				X_start = " << Config.X.start << ", X_end = " << Config.X.end << ", X_Step = " << Config.X.step << ", X_NbStep = " << Config.X.NbStep << "\n";
	file << "				Y_start = " << Config.Y.start << ", Y_end = " << Config.Y.end << ", Y_Step = " << Config.Y.step << ", Y_NbStep = " << Config.Y.NbStep << "\n";
	file << "				Z_start = " << Config.Z.start << ", Z_end = " << Config.Z.end << ", Z_Step = " << Config.Z.step << ", Z_NbStep = " << Config.Z.NbStep << "\n";
	file << "				Root FileOutput = " << Config.nameFile.root << "\n";
	file << "				iterMax = " << Config.Iter.max << "\n";
	file << "				iterMin = " << Config.Iter.min << "\n";
	file << "				rMax = " << Config.rMax << "\n";
	file << "				Filter = " << Config.filter<< "\n";
	file << "				Power = " << Config.power << "\n";
	file << "				dev = " << Config.dev << "\n";
	file << "				IsShow = " << Config.isShow << "\n";
	file << "				NbPoints per step = " << Config.NbPointByStep << "\n";
	file << "				ouput File :  " << Config.nameFile.root << "\n";
	file << "cmd for use this configuration: " << "\n";
	file << "               " << argv[0] << "  ";
	if (Config.W.isFix == 1)
		file << Str_wFix << " " << Config.W.start;
	else
	{
		file << Str_wMin << " " << Config.W.start << " ";
		file << Str_wMax << " " << Config.W.end << " ";
		file << Str_wNbStep << " " << Config.W.NbStep << " ";
	}
	if (Config.X.isFix == 1)
		file << Str_xFix << " " << Config.X.start;
	else
	{
		file << Str_xMin << " " << Config.X.start << " ";
		file << Str_xMax << " " << Config.X.end << " ";
		file << Str_xNbStep << " " << Config.X.NbStep << " ";
	}
	if (Config.Y.isFix == 1)
		file << Str_yFix << " " << Config.Y.start;
	else
	{
		file << Str_yMin << " " << Config.Y.start << " ";
		file << Str_yMax << " " << Config.Y.end << " ";
		file << Str_yNbStep << " " << Config.Y.NbStep << " ";
	}
	if (Config.Z.isFix == 1)
		file << Str_zFix << " " << Config.Z.start;
	else
	{
		file << Str_zMin << " " << Config.Z.start << " ";
		file << Str_zMax << " " << Config.Z.end << " ";
		file << Str_zNbStep << " " << Config.Z.NbStep << " ";
	}
	if (Config.Iter.isFix == 1)
		file << Str_IterFix << " " << Config.Iter.max;
	else
	{
		file << Str_IterMin << " " << Config.Iter.min << " ";
		file << Str_IterMax << " " << Config.Iter.max << " ";
	}
	file << Str_dev << " " << Config.dev << " ";
	file << Str_filter << " " << Config.filter << " ";
	file << Str_power << " " << Config.power << " ";
	file << Str_IsShow << " " << Config.isShow << " ";
	file << Str_rMax << " " << Config.rMax << " ";
	file << Str_Out << " " << Config.nameFile.root << " ";

	file << "\n";
	file.close();

	//Affiche Help
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
	return 0;//Fin du programme
}