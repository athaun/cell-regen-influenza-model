/** 
Viral model with cellular regeneration
Version 2.0

Authors:
Baylor Fain
Asher Haun

Run: nvcc ViralTransmission.cu -o program.out && ./program.out
*/


// Include Header files from same directory
#include "GPU.cuh"
#include "CPU.cuh"
#include "logs.cuh"
#include "simulation.cuh"
#include "interface.cuh"

// Globals to setup the kernals
dim3 BlockConfig, GridConfig;

// Simulation Parameters
int CELL2CELL = 0;
int FREECELL = 1;
float timestep = 0.005;    // Time step for model (No larger than 0.01 hour) 0.005 hr = 18 sec, (1 / 3600) hr = 1 sec
float endtime = 30 * 24; // (2 * 365) * 24;   // in hours
int Save = (1 / timestep); // the number of time the program saves to file, (1 / timestep) results in 1 save every simulated hour
int NumberOfLayers = 607; // 607 is a million hexagons in a circle
int StartRuns = 0;
int NumberOfRuns = 100;

// Physical Parameters
// float MOI[6] = {powf(10,0), powf(10,-1), powf(10,-2), powf(10,-3), powf(10,-4), powf(10,-5)};
float MOI[1] = {powf(10, -2)};
float beta = 2.0; // 2.3 * pow(10,-7); //Infiction rate, units: per hour
float rho = 562800; // 1920
float D = 2.16  * pow(10, -8); // Diffusion rate at 37 degrees celsius unit: m^2 / s //pow(6 * 10,-12) //3.96e - 8
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

// double regenValues[] = {0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0}; // must always contain .0 even if it is an integer value
int totalRegenerations = 0;
double regenParameter;
bool regensAllowed = true;

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
float* timeDead_GPU;
float* th_GPU;
float* ut;
float* ut_GPU;
float* EclipsePhaseLength;
float* EclipsePhaseLength_GPU;
float* InfectionPhaseLength;
float* InfectionPhaseLength_GPU;
float* RegenTime;
float* RegenTime_GPU;
int NumberOfCells;

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
      generates random distribution based off of parameter pasrunSimulation(regenParameter);sed from command line
      (either by individual test run or by multithread.sh)
    */
    // cout << mean << endl;
    random_device rd;
    default_random_engine generator(rd());
    exponential_distribution<double> distribution(mean / 24);
    // cout << (distribution(generator)) << endl;
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
    free(RegenTime);
    free(timeDead);

    if (RUNCPU == 0) {
	    cudaFree(cells_GPU);
	    cudaFree(ecl_GPU);
	    cudaFree(inf_GPU);
	    cudaFree(vtemp_GPU);
	    cudaFree(th_GPU);
	    cudaFree(ut_GPU);
	    cudaFree(EclipsePhaseLength_GPU);
	    cudaFree(InfectionPhaseLength_GPU);
        cudaFree(RegenTime_GPU);
	    cudaFree(state);
        cudaFree(timeDead_GPU);
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

int main(int argc,char *argv[] ) {
    //Checks for Heisenberg status of viral diffusion
    if (D * timestep / pow(deltxprime, 2.0) > 0.5) {
        printf("[FATAL] Change parameters to fit diffusion limits. Value must be under 0.5, current: %.1f. Exiting.", D * timestep / pow(deltxprime, 2.0));
        exit(0);
    }
    //Clear Terminal
    system("clear");

    if(argc > 1) {
       regenParameter = atof(argv[1]); // convert the char to a double
       if (regenParameter != 0) {
         cout << "[SUCCESS] " << regenParameter << " passed from the command line.\n" << endl;
       } else {
         printf("[FATAL] No argument or 0 passed from command line, Exiting.\n");
         exit(0);
       }
     } else {
       printf("[FATAL] No argument passed from command line, Exiting.\n");
     }

    runSimulation(regenParameter);

    printf("[SUCCESS] Program Finished\n");
}
