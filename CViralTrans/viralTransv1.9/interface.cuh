/**

  This interface must contain initializations for all variables and methods that are accessed globally/between files

*/

/** Imports */
#include <stdio.h>
#include <stdlib.h>
#include "cuda.h"
#include <curand_kernel.h>
#include <time.h>
#include <ctime>
#include <math.h>
#include <random>

#include <iomanip>
#include <map>

#include <iostream>
#include <thread>

/** Mehods */
// extern double regenParameter;

void creatingCellLocations();

void runSimulation(float regenParameter);

void deviceSetupAndMemoryAllocation(int Nx, int Ny);

float Te(float TauE, float ne);

float Ti(float TauI, float ni);

float PU1();

float exponentialDistro (double mean);

void errorCheck(const char *message);

void printInitialConditions(int SideLength, int RadiusOfCircle, float regenParameter);

void freeMemory();

using namespace std;

/** variables and constants **/
#define PI 3.1415926535897932f

#define CODETESTINGCONDITIONS 0
#define RUNCPU 0

// Globals to setup the kernals
extern dim3 BlockConfig, GridConfig;

// Simulation Parameters
extern int CELL2CELL;
extern int FREECELL;
extern float timestep;
extern float endtime;
extern int Save;
extern int NumberOfLayers;
extern int StartRuns;
extern int NumberOfRuns;

// extern int SideLength;
// extern int RadiusOfCircle;

// Physical Parameters
// float MOI = pow(10,0); //pow(10,-5) to 1
extern float MOI[1];
extern float beta; // infection rate
extern float rho; // production rate
extern float D; // diffusion coefficient
extern float c; // clearance rate (virus death rate)
extern float deltx; // timestep for diffusion
extern float deltxprime; // part of diffusion "scheme"
extern float Dtsx2; // part of diffusion "scheme"

// Probability Constants
extern float TauI; // avg time cell stays infected
extern float TauE; // avg time cell stays in eclipse stage
extern float ne; // number of eclipse compartments
extern float ni; // number of infected compartments
extern float probi; // Probability per unit time of cell to cell infection (/hour)

extern double regenValues[];
extern int totalRegenerations;
extern double runParameter;
extern double regenParameter;
extern bool regensAllowed;

//Global Variables
extern char Path_to_Folder[100];
extern char directory[100];
extern char** LocationData;
extern char* cells;
extern char* cells_GPU;
extern float* ecl;
extern float* ecl_GPU;
extern float* inf;
extern float* inf_GPU;
extern float* vtemp;
extern float* vtemp_GPU;
extern float* th;
extern float* timeDead;
extern float* timeDead_GPU;
extern float* th_GPU;
extern float* ut;
extern float* ut_GPU;
extern float* EclipsePhaseLength;
extern float* EclipsePhaseLength_GPU;
extern float* InfectionPhaseLength;
extern float* InfectionPhaseLength_GPU;
extern float* RegenTime;
extern float* RegenTime_GPU;
extern int NumberOfCells;
extern int NumberDead;

extern int NumberDead1;
extern int NumberInfected1;
extern int NumberEclipse1;
extern int NumberHealthy1;
extern float AmountOfVirus;

extern curandState *state;

extern int NumberHealthy;
extern int NumberEclipse;
extern int NumberInfected;
extern int NumberDead;
extern int NumberVirus;
