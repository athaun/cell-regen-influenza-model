
//nvcc ViralTransmission.cu -o program.out && ./program.out

/**
  NOTES:
  - figure out multithreaded aproach for the values of regeneration
  - maybe do ^^^ for MOI???

  CHANGELOG:
  1.12.2020
  - separated ViralTransmission.cu into multiple files, and added an interface so that they can communicate
  - removed the old cerialViralTransmission() method
  - started to rename variables (especially in for loops) to names that make more sense, ie from j to y when sorting through a 2D array
  1.13.2020
  - finished the interface with method declarations
  - cleaned up if/else ladders with switch statements
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
#include <iostream>

// Include Header files from same directory
#include "GPU.cuh"
#include "CPU.cuh"
#include "logs.cuh"
#include "interface.cuh"

using namespace std;

// Globals to setup the kernals
dim3 BlockConfig, GridConfig;

// Simulation Parameters
int CELL2CELL = 1;
int FREECELL = 0;
float timestep = 0.005;    // Time step for model (No larger than 0.01 hour) 0.005 hr = 18 sec, (1 / 3600) hr = 1 sec
float endtime = (2 * 365) * 24;   // Days in hours
int Save = (1 / timestep); // the number of time the program saves to file, (1 / timestep) results in 1 save every simulated hour
int NumberOfLayers = 7; // 607 is a million hexagons in a circle
int StartRuns = 0;
int NumberOfRuns = 1;

// Physical Parameters
// float MOI = pow(10,0); //pow(10,-5) to 1
float beta = 2.0; // 2.3 * pow(10,-7); //Infiction rate, units: per hour
float rho = 562800; // 1920
float D = 6 * pow(10, -12); // Diffusion rate at 37 degrees celsius unit: m^2 / s //pow(6 * 10,-12) //3.96e - 8
float c = 0.105; // Clearance rate, units: per hour
float deltx = 25.0 * pow(10, -6);
float deltxprime = deltx * 2;
float Dtsx2 = D * timestep * pow(deltxprime, -2);

// Probability Constants
float TauI = 12.0;  // Avg time for infection
float TauE = 6.0;   // Avg time for eclipse
float ne = 30.0;    // Number of eclipse compartments?
float ni = 100.0;   // Number of infected compartments?
float probi = 0.2;  // Probability per unit time of cell to cell infection (/hour)

double regenValues[] = {0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0}; // must always contain .0 even if it is an integer value
int totalRegenerations = 0;
double runParameter;

//Global Variables
char Path_to_Folder[100] = "";
char directory[100] = "";
char** LocationData;
char* cells;
char* cells_GPU;
float* ecl;
float* ecl_GPU;
float* inf;
float* inf_GPU;
float* vtemp;
float* vtemp_GPU;
float* th;
float* timeDead;
float* th_GPU;
float* ut;
float* ut_GPU;
float* EclipsePhaseLength;
float* EclipsePhaseLength_GPU;
float* InfectionPhaseLength;
float* InfectionPhaseLength_GPU;
int NumberOfCells;
int NumberDead;

int NumberDead1;
int NumberInfected1;
int NumberEclipse1;
int NumberHealthy1;
float AmountOfVirus;

curandState *state;

/** Functions */

float Te(float TauE, float ne) {
    /**
    Picks a random number from the gamma distribution
    The number is to be used as a time step in the Eclipse Time Matrix
    */
    random_device rd;
    default_random_engine generator(rd());
    gamma_distribution<double> distribution(TauE, TauE / sqrt(ne));

    return distribution(generator);
}

float Ti(float TauI, float ni) {
    /**
    Picks a random number from the gamma distribution
    The number is to be used as a time step in the Infected Time Matrix
    */
    random_device rd;
    default_random_engine generator(rd());
    gamma_distribution<double> distribution(TauI, TauI / sqrt(ni));

    return distribution(generator);
}

float PU1() {
    // Picks a random number from a uniform distribution
    random_device rd;
    default_random_engine generator(rd());
    uniform_real_distribution<double> distribution(0.0,1.0);

    return distribution(generator);
}

float exponentialDistro (double mean) {

    /**
      generates random distribution based off of parameter passed from command line
      (either by individual test run or by multithread.sh)
    */
    runParameter = mean;

    random_device rd;
    default_random_engine generator(rd());
    exponential_distribution<double> distribution(runParameter * 24);

    return distribution(generator);
}

void freeMemory() {
    for (int i = 0; i < ((2 * NumberOfLayers)-1); i ++) {
        free(LocationData[i]);
    }
    free(LocationData);
    free(cells);
    free(ecl);
    free(inf);
    free(vtemp);
    free(th);
    free(ut);
    free(EclipsePhaseLength);
    free(InfectionPhaseLength);

    if (RUNCPU == 0) {
	    cudaFree(cells_GPU);
	    cudaFree(ecl_GPU);
	    cudaFree(inf_GPU);
	    cudaFree(vtemp_GPU);
	    cudaFree(th_GPU);
	    cudaFree(ut_GPU);
	    cudaFree(EclipsePhaseLength_GPU);
	    cudaFree(InfectionPhaseLength_GPU);
	    cudaFree(state);
    }
}

void errorCheck(const char *message) {
  cudaError_t error;
  error = cudaGetLastError();

  if (error != cudaSuccess) {
    printf("\n CUDA ERROR: %s = %s\n", message, cudaGetErrorString(error));
    exit(0);
  }
}

int main(void) {
    //Checks for Heisenberg status of viral diffusion
    if (D * timestep / pow(deltxprime, 2.0) > 0.5) {
        printf("%.1f", D * timestep / pow(deltxprime, 2.0));
        printf("CHANGE PARAMETERS TO FIT DIFFUSION LIMITS. VALUE MUST BE UNDER 0.5. VALUE SHOWN ABOVE");
        exit(0);
    }
    //Clear Terminal
    system("clear");

    // float MOI[6] = {powf(10,0), powf(10,-1), powf(10,-2), powf(10,-3), powf(10,-4), powf(10,-5)};
    float MOI[1] = {powf(10, -2)};
    for(int a = 0; a < sizeof(regenValues)/sizeof(regenValues[0]); a ++) {
        for (int q = 0; q < (sizeof(MOI) / sizeof(MOI[0])); q ++) {
            for (int BigIndex = 0; BigIndex < NumberOfRuns; BigIndex ++) {
              // run the test as many times as necessary to get through all of the regen parameters, MOI values, and runs defined in NumberOfRuns

                printf("\nStarting run %d\n", (BigIndex + 1));

                //Creating Save Path
                creatingPathToFolderAndDirectory(StartRuns + BigIndex, NumberOfLayers, MOI[q]);
                //Creating placeholder variables for multiple runs
                int cell2cell = CELL2CELL;
                int freecell = FREECELL;

                //Building Cells
                creatingCellLocations();

                //Number of Cells
                //Number of initial infected cells
                int Ni = NumberOfCells * MOI[q]; if (Ni < 1) {
                  printf("Use larger MOI\n");
                  exit(0);
                }
                int Nx = (2 * NumberOfLayers - 1);      // Range of cells on x axis
                int Ny = (2 * NumberOfLayers - 1);      // Range of cells on y axis

                //Making empty matrices
                allocateMemory(Nx, Ny);

                //Initializing
                initailConditions(Nx, Ny);

                //Deletes files and initial with values
                if (BigIndex == 0) {
                    printToFileCellAndVirusInitial(Nx, Ny, NumberOfLayers);
                }

                printToFileCellAndVirusAnalysisInitial(Nx, Ny);

                //Infects a random cell, now seen as (e)
                infectRandomCells(Nx, Ny, Ni);

                if (RUNCPU == 0) {
                    // Run GPU code
                    cudaMalloc((void**)&state, Nx * Ny * sizeof(int));
                    errorCheck("cudaMalloc Random Setup");
                    cuRand_Setup<<<GridConfig, BlockConfig>>>(state);
                    errorCheck("Random Setup");

                    loadConstants(MOI[q]);

                    deviceSetupAndMemoryAllocation(Nx, Ny);

          	        cudaMemcpy(cells_GPU, cells, Nx * Ny * 2 * sizeof(char), cudaMemcpyHostToDevice );
          	        errorCheck("cudaMemcpy cells HtoD");
          	        cudaMemcpy(vtemp_GPU, vtemp, Nx * Ny * 2 * sizeof(float), cudaMemcpyHostToDevice );
          	        errorCheck("cudaMemcpy vtemp HtoD");
          	        cudaMemcpy(ut_GPU, ut, Nx * Ny * sizeof(float), cudaMemcpyHostToDevice );
          	        errorCheck("cudaMemcpy ut HtoD");

          	        cudaMemcpy(ecl_GPU, ecl, Nx * Ny * sizeof(float), cudaMemcpyHostToDevice );
          	        errorCheck("cudaMemcpy ecl HtoD");
          	        cudaMemcpy(inf_GPU, inf, Nx * Ny * sizeof(float), cudaMemcpyHostToDevice );
          	        errorCheck("cudaMemcpy inf HtoD");
          	        cudaMemcpy(th_GPU, th, Nx * Ny * sizeof(float), cudaMemcpyHostToDevice );
          	        errorCheck("cudaMemcpy th HtoD");

          	        cudaMemcpy(EclipsePhaseLength_GPU, EclipsePhaseLength, Nx * Ny * sizeof(float), cudaMemcpyHostToDevice );
          	        errorCheck("cudaMemcpy EclipsePhaseLength HtoD");
          	        cudaMemcpy(InfectionPhaseLength_GPU, InfectionPhaseLength, Nx * Ny * sizeof(float), cudaMemcpyHostToDevice );
          	        errorCheck("cudaMemcpy InfectionPhaseLength HtoD");
                }

                //Runs simulation
                int NumberofTimeSteps = endtime / timestep;
                int NumberofSavedTimeSteps = NumberofTimeSteps / Save;
                int timestepcount = 0;    //equal to the number of ts elapsed
                while(timestepcount < (NumberofTimeSteps - 1)) {

                    if (RUNCPU == 0) {
                        // GPU code
                        kernel<<<GridConfig,BlockConfig>>>(cells_GPU, vtemp_GPU, ut_GPU, ecl_GPU, inf_GPU, th_GPU, EclipsePhaseLength_GPU, InfectionPhaseLength_GPU, SystemConstants, cell2cell, freecell, state, NumberOfLayers);
                    } else {
                        modifiedCerialViralTransmission(Nx, Ny, cell2cell, freecell);
                    }

                    if ((timestepcount % Save) == 0) {
                        if (RUNCPU == 0) {
                            cudaMemcpy(cells, cells_GPU, Nx * Ny * 2 * sizeof(char), cudaMemcpyDeviceToHost );
                            errorCheck("cudaMemcpy cells DtoH");
                            cudaMemcpy(vtemp, vtemp_GPU, Nx * Ny * 2 * sizeof(float), cudaMemcpyDeviceToHost );
                            errorCheck("cudaMemcpy vtemp DtoH");
                        }

                        // Analysis dish
                        NumberDead1 = 0;
                        NumberInfected1 = 0;
                        NumberEclipse1 = 0;
                        NumberHealthy1 = 0;
                        AmountOfVirus = 0.0;
                        for (int y = 0; y < Ny; y ++) {
                            for (int x = 0; x < Nx; x ++) {
                                AmountOfVirus += vtemp[x + (Nx * y) + (Nx * Ny * 0)];
                                switch (cells[x + Nx * y + Nx * Ny * 0]) { // why multiply by 0?
                                  case 'h':
                                    NumberHealthy1 ++;
                                    break;
                                  case 'e':
                                    NumberEclipse1 ++;
                                    break;
                                  case 'i':
                                    NumberInfected1 ++;
                                    break;
                                  case 'd':
                                    NumberDead1 ++;
                                    break;
                                }
                            }
                        }

                        // Prints status of cells virus
                        if (BigIndex == 0) {
                            printToFileCellAndVirus(Nx, Ny, NumberOfLayers);
                        }

                        printToFileCellAndVirusAnalysis(timestepcount * timestep);
                    }

                    //Number of days completed
                    if ((timestepcount % (24 * int(1 / timestep))) == 0) {
                        printf("%.0f Day\n", (timestepcount * timestep) / 24);
                    }

                    if ((NumberHealthy1 == 0)) {
                        cell2cell = 0;
                        freecell = 0;
                    } else {
                        cell2cell = CELL2CELL;
                        freecell = FREECELL;
                    }

                    //End Code if Virus has below 10
                    if ((AmountOfVirus < pow(10, 1.0)) && (NumberDead1 == NumberOfCells)) {
                        break;
                    }

                    timestepcount ++;
                }

                //Writes a file with all of our parameters / variables
                createParameterFile(timestep, NumberofSavedTimeSteps, endtime, timestepcount, AmountOfVirus, rho, D, deltxprime, c, probi);

                printf("\n%d of %d Runs Done\n", (BigIndex + 1), NumberOfRuns);

                freeMemory();
            }
        }
    }
    printf("PROGRAM DONE\n");
}
