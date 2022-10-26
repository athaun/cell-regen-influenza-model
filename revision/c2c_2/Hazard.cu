//##############################################################################
//#                                                                            #
//#                          Virus Model                                       #
//#                                                                            #
//##############################################################################

// nvcc Hazard_double.cu -o program.out && ./program.out

// Number of Layer, Number of Cells
// 7              , 121
// 20             , 1027
// 61             , 10009
// 193            , 100969
// 430            , 501535
// 607            , 1001365

/*
    Using Camal case:
        C Functions start with lower case
        Variables start with upper case
*/

#include <curand_kernel.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <ctime>
#include <iostream>
#include <random>

#include "cuda.h"
#include <chrono>

using namespace std;

#define PI 3.1415926535897932f

#define CODETESTINGCONDITIONS 0
#define RUNCPU 0

// Globals to setup the kernals
dim3 BlockConfig, GridConfig;

// Simulation Parametersf
int CELL2CELL = 1;
int FREECELL = 0;

// Change initial conditions
int INITIALCELL = 1;
int INITIALVIRUS = 0;

// 0.000625, 0.00128, 0.0025, 0.005, 0.01, 0.02, 0.04
float timestep = 0.005;     // Time step for model (No larger than 0.01 hour) 0.005 hr = 18 sec, (1/3600) hr = 1 sec
float endtime = 100*24;   // Days in hours 280;//
int Save = (1 / timestep);  // the number of time the program saves to file, (1/timestep) results in 1 save every simulated hour
int NumberOfLayers = 607;   // 7 //20 //61 //193 //607 is a million hexigon in a circle
int StartRuns = 0;
int NumberOfRuns = 100;

// Physical Parameters
// float MOI = pow(10,0); //pow(10,-5) to 1
float _beta = 27.0;              // 2.3*pow(10,-7); //Infiction rate, units: per hour
float rho = 3000.0;           // 562800
float D = 2.16 * pow(10, -8);  // Diffusion rate at 37 degrees celsius unit: m^2/h
float c = 0.25;                 // 0.105 Clearance rate, units: per hour
float deltx = 25.0 * pow(10, -6);
float deltxprime = deltx * 2;
float Dtsx2 = D * timestep * pow(deltxprime, -2);

// Probability Constants
float TauI = 26.0;  // Avg time for infection
float TauE = 16.0;   // Avg time for eclipse
float ne = 30.0;    // Number of eclipse compartments?
float ni = 100.0;   // Number of infected compartments?
// float probi = 0.2;  //Probability per unit time of cell to cell infection (/hour)

// double regenValues[] = {0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0}; // must always contain .0 even if it is an integer value
int totalRegenerations = 0;
double regenParameter;
bool regensAllowed = true;

// Global Variables
char Path_to_Folder[100] = "";
char directory[100] = "";
char **LocationData;
char *cells;
char *cells_GPU;
float *ecl;
float *ecl_GPU;
float *inf;
float *inf_GPU;
double *vtemp;
double *vtemp_GPU;
float *th;
float *timeDead;
float *timeDead_GPU;
float *th_GPU;
float *ut;
float *ut_GPU;
float *EclipsePhaseLength;
float *EclipsePhaseLength_GPU;
float *InfectionPhaseLength;
float *InfectionPhaseLength_GPU;
float *RegenTime;
float *RegenTime_GPU;
int NumberOfCells;

int NumberDead1;
int NumberInfected1;
int NumberEclipse1;
int NumberHealthy1;
float AmountOfVirus;

curandState *state;

// Functions
float Te(float TauE, float ne) {
    //    Picks a random number from the gamma distribution
    //    The number is to be used as a time step in the Eclipse Time Matrix
    random_device rd;
    default_random_engine generator(rd());
    gamma_distribution<double> distribution(ne, TauE / ne);

    return distribution(generator);
}

float Ti(float TauI, float ni) {
    //    Picks a random number from the gamma distribution
    //    The number is to be used as a time step in the Infected Time Matrix
    random_device rd;
    default_random_engine generator(rd());
    gamma_distribution<double> distribution(ni, TauI / ni);

    return distribution(generator);
}

float PU1() {
    //    Picks a random number from a uniform distribution
    //    This probability
    random_device rd;
    default_random_engine generator(rd());
    uniform_real_distribution<double> distribution(0.0, 1.0);

    return distribution(generator);
}

float exponentialDistro(double mean) {
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

void creatingPathToFolderAndDirectory(int BigIndex, int NumberOfLayers, float MOI, float regenParameter) {
    char TransmissionType[10] = "";
    if (CELL2CELL == 1) {
        if (FREECELL == 1) {
            strcat(TransmissionType, "Both");
        } else {
            strcat(TransmissionType, "CELL2CELL");
        }
    } else if (CELL2CELL == 0) {
        if (FREECELL == 0) {
            strcat(TransmissionType, "Neither");
        } else {
            strcat(TransmissionType, "FREECELL");
        }
    }

    char Buffer[5];  // Buffer String For Conversion To Char
    char TheCurrentTime[50];
    time_t RawTime = time(NULL);
    tm *SpecificMoment = localtime(&RawTime);

    strcpy(Path_to_Folder, "");
    strcpy(directory, "");

    if (RUNCPU == 1) {
        strcat(Path_to_Folder, "ViralModel/");
    } else {
        strcat(Path_to_Folder, "ViralModel/");
    }
    // strftime(TheCurrentTime, 50, "%m-%d/%I:%M", SpecificMoment);
    strftime(TheCurrentTime, 50, "%m-%d/", SpecificMoment);
    strcat(Path_to_Folder, TheCurrentTime);
    // strcat(Path_to_Folder,"_");
    sprintf(Buffer, "%d", NumberOfLayers);
    strcat(Path_to_Folder, Buffer);
    strcat(Path_to_Folder, "_");
    sprintf(Buffer, "%d", BigIndex);
    strcat(Path_to_Folder, Buffer);
    strcat(Path_to_Folder, "-");
    strcat(Path_to_Folder, TransmissionType);
    strcat(Path_to_Folder, "_");

    sprintf(Buffer, "%.1f", log10(MOI));
    strcat(Path_to_Folder, Buffer);
    strcat(Path_to_Folder, "-");
    strcat(Path_to_Folder, "MOI");

    strcat(Path_to_Folder, "_");
    sprintf(Buffer, "%.3f", regenParameter);
    strcat(Path_to_Folder, Buffer);
    strcat(Path_to_Folder, "-RP");
    // regenParameter

    strcat(directory, "mkdir -p ");
    strcat(directory, Path_to_Folder);
    int check = system(strdup(directory));
    if (check != 0) {
        exit(0);
    }
}

void creatingCellLocations() {
    float SideLenght = (2.0 / 3.0);
    int RadiusScale = 0;
    for (int i = 0; i < NumberOfLayers; i++) {
        if (i == 0) {
            RadiusScale = RadiusScale + 1;
        } else {
            if ((i) % 2 == 1) {
                RadiusScale = RadiusScale + 1;
            } else {
                RadiusScale = RadiusScale + 2;
            }
        }
    }
    float RadiusOfCircle = SideLenght * RadiusScale;

    int count = 0;
    for (int i = 0; i < NumberOfLayers; i++) {
        count = count + i;
    }
    int NumberOfHexagons = (count)*6 + 1;

    float **coord;
    int n = NumberOfHexagons;
    int m = 3;
    coord = (float **)calloc(n, sizeof(float *));
    for (int i = 0; i < n; i++) {
        coord[i] = (float *)calloc(m, sizeof(float));
    }

    float **percyclecoord;
    n = NumberOfHexagons;
    m = 3;
    percyclecoord = (float **)calloc(n, sizeof(float *));
    for (int i = 0; i < n; i++) {
        percyclecoord[i] = (float *)calloc(m, sizeof(float));
    }

    int temp;
    for (int j = 0; j < NumberOfLayers; j++) {
        for (int i = 0; i < (2 * j); i++) {
            if (i < j) {
                temp = i;
            }
            percyclecoord[i + (j - 1) * j + 1][0] = -temp - 1;
            percyclecoord[i + (j - 1) * j + 1][1] = temp + j - i;
            percyclecoord[i + (j - 1) * j + 1][2] = -j + 1 + i;
        }
    }
    float c0[3] = {percyclecoord[0][0], percyclecoord[0][1], percyclecoord[0][2]};
    coord[0][2] = c0[2];
    coord[0][1] = c0[1];
    coord[0][0] = c0[0];

    count = 0;
    for (int j = 0; j < (NumberOfHexagons / 3); j++) {
        for (int i = 0; i < 3; i++) {
            coord[(i + 0) % 3 + 3 * j + 1][2] = percyclecoord[j + 1][i] + c0[i];
            coord[(i + 1) % 3 + 3 * j + 1][1] = percyclecoord[j + 1][i] + c0[i];
            coord[(i + 2) % 3 + 3 * j + 1][0] = percyclecoord[j + 1][i] + c0[i];
        }
    }

    float hi = coord[0][0];
    float vi = coord[0][2];
    float xmin = INFINITY;
    float xcoord;
    float ycoord;
    double dist;
    for (int i = 0; i < NumberOfHexagons; i++) {
        xcoord = coord[i][0];
        if (coord[i][0] < xmin) {
            xmin = coord[i][0];
        }
        ycoord = (2.0 * sin(PI * (60.0 / 180.0)) * (coord[i][1] - coord[i][2]) / 3.0) + vi;
        dist = sqrtf(pow(double(xcoord - hi), 2.0) + pow(double(ycoord - vi), 2.0));
        if (dist >= RadiusOfCircle) {
            coord[i][0] = 5000.0;
            coord[i][1] = 0.0;
            coord[i][2] = 0.0;
        }
    }

    n = ((2 * NumberOfLayers) - 1);
    m = ((2 * NumberOfLayers) - 1);
    LocationData = (char **)malloc(n * sizeof(char *));
    for (int j = 0; j < n; j++) {
        LocationData[j] = (char *)malloc(m * sizeof(char));
        for (int i = 0; i < m; i++) {
            LocationData[j][i] = 'o';
        }
    }

    NumberOfCells = 0;
    for (int i = 0; i < NumberOfHexagons; i++) {
        if (coord[i][0] != 5000.0) {
            LocationData[int(coord[i][2]) - int(xmin)][int(coord[i][0]) - int(xmin)] = 'h';
            NumberOfCells = NumberOfCells + 1;
        }
    }

    //    char File1[100] = "";
    //    strcat(File1,Path_to_Folder);
    //    strcat(File1,"/InitialCellLocations.txt");
    //    FILE *outfile1 = fopen(File1,"a");
    //    if (outfile1 == NULL){
    //        printf("Error opening file!\n");
    //        exit(0);
    //    }
    //
    //    for(int i=0; i<((2*NumberOfLayers)-1); i++){
    //        for(int j=0; j<((2*NumberOfLayers)-1); j++){
    //            fprintf(outfile1,"%c,",LocationData[i][j]);
    //        }
    //            fprintf(outfile1,"\n");
    //    }
    //    fclose(outfile1);

    char File2[100] = "";
    strcat(File2, Path_to_Folder);
    strcat(File2, "/Parameters.txt");
    FILE *outfile2 = fopen(File2, "w");
    if (outfile2 == NULL) {
        printf("Error opening file!\n");
        exit(0);
    }
    fprintf(outfile2, "Hexagon Side Lenght = %f\n", SideLenght);
    fprintf(outfile2, "Number of Layers = %d\n", NumberOfLayers);
    fprintf(outfile2, "Radius of Circle = %f\n", RadiusOfCircle);
    fprintf(outfile2, "Number of Cells = %d\n", NumberOfCells);
    fprintf(outfile2, "Regeneration rate = %d\n", regenParameter);
    fclose(outfile2);

    for (int i = 0; i < NumberOfHexagons; i++) {
        free(coord[i]);
    }
    free(coord);

    for (int i = 0; i < NumberOfHexagons; i++) {
        free(percyclecoord[i]);
    }
    free(percyclecoord);
}

void allocateMemory(int Nx, int Ny) {
    // Produces a matrix for the cells
    cells = (char *)malloc(Nx * Ny * 2 * sizeof(char));
    // Produces a matrix that will track the amount virus above each cell
    vtemp = (double *)malloc(Nx * Ny * 2 * sizeof(double));
    // Produces a univeral time matrix (ut)
    ut = (float *)malloc(Nx * Ny * sizeof(float));

    // Produces a time matrix for after eclipse phase (e)
    ecl = (float *)malloc(Nx * Ny * sizeof(float));
    // Produces a time matrix for after infection phase (i)
    inf = (float *)malloc(Nx * Ny * sizeof(float));
    // Produces a time matrix hor healthy cells (t)
    th = (float *)malloc(Nx * Ny * sizeof(float));

    timeDead = (float*) malloc(Nx * Ny * sizeof(float));

    // Produces an array of eclipse phase durations for cells
    EclipsePhaseLength = (float *)malloc(Nx * Ny * sizeof(float));
    // Produces an array of infection phase durations for cells
    InfectionPhaseLength = (float *)malloc(Nx * Ny * sizeof(float));

    RegenTime = (float*) malloc(Nx * Ny * sizeof(float));
}

void initailConditions(int Nx, int Ny, double MOI) {
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            for (int k = 0; k < 2; k++) {
                cells[i + Nx * j + Nx * Ny * k] = LocationData[i][j];
                vtemp[i + Nx * j + Nx * Ny * k] = 0.0;
            }
            ut[i + Nx * j] = 0.0;
            ecl[i + Nx * j] = 0.0;
            inf[i + Nx * j] = 0.0;
            th[i + Nx * j] = 0.0;
            timeDead[i + Nx * j] = 0.0;
            EclipsePhaseLength[i + Nx * j] = Te(TauE, ne);
            InfectionPhaseLength[i + Nx * j] = Ti(TauI, ni);
            RegenTime[i + Nx * j] = exponentialDistro(regenParameter);
        }
    }
    if (INITIALVIRUS == 1) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                for (int k = 0; k < 2; k++) {
                    if (cells[i + Nx * j + Nx * Ny * k] != 'o') {
                        vtemp[i + Nx * j + Nx * Ny * k] = MOI;
                    }
                }
            }
        }
    }
}

void infectANumberOfCellsRandomly(int Nx, int Ny, int Ni) {
    if (CODETESTINGCONDITIONS == 1) {
        cells[(NumberOfLayers - 1) + Nx * (NumberOfLayers - 1) + Nx * Ny * 0] = 'i';
        cells[(NumberOfLayers - 1) + Nx * (NumberOfLayers - 1) + Nx * Ny * 1] = 'i';  // Only the center cell
    } else {
        srand(time(NULL));
        int randx;
        int randy;
        int NumberOfInfectedCellsCount = 0;
        while (NumberOfInfectedCellsCount < Ni) {
            randx = (rand() % Nx);
            randy = (rand() % Ny);
            if ((cells[randx + Nx * randy + Nx * Ny * 0] != 'o') && (cells[randx + Nx * randy + Nx * Ny * 0] == 'h')) {
                cells[randx + Nx * randy + Nx * Ny * 0] = 'e';
                cells[randx + Nx * randy + Nx * Ny * 1] = 'e';
                ecl[randx + Nx * randy] = Te(TauE, ne);
                NumberOfInfectedCellsCount ++;
            }
        }
    }
}

void printToFileCellAndVirusInitial(int Nx, int Ny, int NumberOfLayers) {
    char File3[100] = "";
    strcat(File3, Path_to_Folder);
    strcat(File3, "/cells_over_time.txt");
    FILE *outfile3 = fopen(File3, "w");
    if (outfile3 == NULL) {
        printf("Error opening file!\n");
        exit(0);
    }
    for (int i = 0; i < ((2 * NumberOfLayers) - 1); i++) {
        for (int j = 0; j < ((2 * NumberOfLayers) - 1); j++) {
            fprintf(outfile3, "%c,", LocationData[i][j]);
        }
        fprintf(outfile3, "\n");
    }
    fclose(outfile3);

    char File4[100] = "";
    strcat(File4, Path_to_Folder);
    strcat(File4, "/virus_over_time.txt");
    FILE *outfile4 = fopen(File4, "w");
    if (outfile4 == NULL) {
        printf("Error opening file!\n");
        exit(0);
    }
    for (int i = 0; i < ((2 * NumberOfLayers) - 1); i++) {
        for (int j = 0; j < ((2 * NumberOfLayers) - 1); j++) {
            fprintf(outfile4, "%f,", 0.0);
        }
        fprintf(outfile4, "\n");
    }
    fclose(outfile4);
}

void printToFileCellAndVirusAnalysisInitial (int Nx, int Ny) {
    NumberDead1 = 0;
    NumberInfected1 = 0;
    NumberEclipse1 = 0;
    NumberHealthy1 = 0;
    AmountOfVirus = 0.0;
    for (int j = 0; j < Ny; j ++) {
        for (int i = 0; i < Nx; i ++) {
            AmountOfVirus = AmountOfVirus + vtemp[i + Nx * j + Nx * Ny * 0];

            if (cells[i + Nx * j + Nx * Ny * 0] == 'd') {
                NumberDead1 ++;
            }
            else if (cells[i + Nx * j + Nx * Ny * 0] == 'i') {
                NumberInfected1 ++;
            }
            else if (cells[i + Nx * j + Nx * Ny * 0] == 'e') {
                NumberEclipse1 ++;
            }
            else if (cells[i + Nx * j + Nx * Ny * 0] == 'h') {
                NumberHealthy1 ++;
            }
        }
    }

    char File9[100] = "";
    strcat(File9,Path_to_Folder);
    strcat(File9,"/PerTimeStep.txt");
    FILE *outfile9 = fopen(File9,"w");
    if (outfile9 == NULL) {
        printf("Error opening file!\n");
        exit(0);
    }

    fprintf(outfile9,"%0.0f, %d, %d, %d, %d, %f, %d", 0.0, NumberHealthy1, NumberEclipse1, NumberInfected1, NumberDead1, AmountOfVirus, totalRegenerations);
    fprintf(outfile9,"\n");

    fclose(outfile9);
}

void printToFileCellAndVirus(int Nx, int Ny, int NumberOfLayers) {
    char File5[100] = "";
    strcat(File5, Path_to_Folder);
    strcat(File5, "/cells_over_time.txt");
    FILE *outfile5 = fopen(File5, "a");
    if (outfile5 == NULL) {
        printf("Error opening file!\n");
        exit(0);
    }
    for (int i = 0; i < ((2 * NumberOfLayers) - 1); i++) {
        for (int j = 0; j < ((2 * NumberOfLayers) - 1); j++) {
            fprintf(outfile5, "%c,", cells[i + Nx * j + Nx * Ny * 0]);
        }
        fprintf(outfile5, "\n");
    }
    fclose(outfile5);

    char File6[100] = "";
    strcat(File6, Path_to_Folder);
    strcat(File6, "/virus_over_time.txt");
    FILE *outfile6 = fopen(File6, "a");
    if (outfile6 == NULL) {
        printf("Error opening file!\n");
        exit(0);
    }
    for (int i = 0; i < ((2 * NumberOfLayers) - 1); i++) {
        for (int j = 0; j < ((2 * NumberOfLayers) - 1); j++) {
            fprintf(outfile6, "%f,", vtemp[i + Nx * j + Nx * Ny * 1]);
        }
        fprintf(outfile6, "\n");
    }
    fclose(outfile6);
}

void printToFileCellAndVirusAnalysis(float timestep) {
    char File8[100] = "";
    strcat(File8, Path_to_Folder);
    strcat(File8, "/PerTimeStep.txt");
    FILE *outfile8 = fopen(File8, "a");
    if (outfile8 == NULL) {
        printf("Error opening file!\n");
        exit(0);
    }

    fprintf(outfile8,"%0.0f, %d, %d, %d, %d, %f, %d", timestep + 1, NumberHealthy1, NumberEclipse1, NumberInfected1, NumberDead1, AmountOfVirus, totalRegenerations);
    fprintf(outfile8, "\n");

    fclose(outfile8);
}

void createParameterFile(float timestep, int NumberofSavedTimeSteps, float endtime, float timestepcount, double AmountOfVirus, float rho, float D, float deltxprime, float c, float probi) {
    char File7[100] = "";
    strcat(File7, Path_to_Folder);
    strcat(File7, "/Parameters.txt");
    FILE *outfile7 = fopen(File7, "a");
    if (outfile7 == NULL) {
        printf("Error opening file!\n");
        exit(0);
    }
    fprintf(outfile7, "Time Step = %f\n", timestep);
    fprintf(outfile7, "Number of Saved Time Steps = %d\n", NumberofSavedTimeSteps);
    fprintf(outfile7, "Initial End Time = %f\n", endtime);
    fprintf(outfile7, "Actual Hours Simulated = %f\n", timestepcount * timestep);
    fprintf(outfile7, "Final Amount of Virus = %f\n", AmountOfVirus);
    fprintf(outfile7, "rho = %f\n", rho);
    fprintf(outfile7, "log10(D) = %f\n", log10(D));
    fprintf(outfile7, "delta x = %f\n", deltxprime);
    fprintf(outfile7, "c = %f\n", c);
    fprintf(outfile7, "Probability of cell to cell infection: %f\n", probi);
    fclose(outfile7);
}

void freeMemory() {
    for (int i = 0; i < ((2 * NumberOfLayers) - 1); i++) {
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

struct systemConstantsStruct {
    double MOI;
    float beta;
    float rho;
    float D;
    float c;
    float deltx;
    float deltxprime;
    float Dtsx2;

    float TauI;
    float TauE;
    float ne;
    float ni;
    float probi;

    float timestep;
};

systemConstantsStruct SystemConstants;

void loadConstants(double MOI, float probi) {
    SystemConstants.MOI = MOI;
    SystemConstants.beta = _beta;
    SystemConstants.rho = rho;
    SystemConstants.D = D;
    SystemConstants.c = c;
    SystemConstants.deltx = deltx;
    SystemConstants.deltxprime = deltxprime;
    SystemConstants.Dtsx2 = Dtsx2;

    SystemConstants.TauI = TauI;
    SystemConstants.TauE = TauE;
    SystemConstants.ne = ne;
    SystemConstants.ni = ni;
    SystemConstants.probi = probi;

    SystemConstants.timestep = timestep;
}

void deviceSetupAndMemoryAllocation(int Nx, int Ny) {
    BlockConfig.x = 16;
    BlockConfig.y = 16;
    BlockConfig.z = 1;

    GridConfig.x = (Nx - 1) / BlockConfig.x + 1;
    GridConfig.y = (Ny - 1) / BlockConfig.y + 1;
    GridConfig.z = 1;

    cudaMalloc((void **)&cells_GPU, Nx * Ny * 2 * sizeof(char));
    errorCheck("cudaMalloc cells Mem");
    cudaMalloc((void **)&vtemp_GPU, Nx * Ny * 2 * sizeof(double));
    errorCheck("cudaMalloc vtemp Mem");
    cudaMalloc((void **)&ut_GPU, Nx * Ny * sizeof(float));
    errorCheck("cudaMalloc ut Mem");

    cudaMalloc((void **)&ecl_GPU, Nx * Ny * sizeof(float));
    errorCheck("cudaMalloc ecl Mem");
    cudaMalloc((void **)&inf_GPU, Nx * Ny * sizeof(float));
    errorCheck("cudaMalloc inf Mem");
    cudaMalloc((void **)&th_GPU, Nx * Ny * sizeof(float));
    errorCheck("cudaMalloc th Mem");

    cudaMalloc((void **)&EclipsePhaseLength_GPU, Nx * Ny * sizeof(float));
    errorCheck("cudaMalloc EclipsePhaseLength Mem");
    cudaMalloc((void **)&InfectionPhaseLength_GPU, Nx * Ny * sizeof(float));
    errorCheck("cudaMalloc InfectionPhaseLength Mem");

    cudaMalloc((void**)&RegenTime_GPU, Nx * Ny * sizeof(float));
	errorCheck("cudaMalloc RegenTime Mem");
    cudaMalloc((void**)&timeDead_GPU, Nx * Ny * sizeof(float));
	errorCheck("cudaMalloc timeDead Mem");
}

__global__ void cuRand_Setup(curandState *state) {
    int Row = threadIdx.x + blockIdx.x * blockDim.x;
    int Column = threadIdx.y + blockIdx.y * blockDim.y;
    int offsetx = blockDim.x * gridDim.x;

    int id = Row + offsetx * Column;
    curand_init(clock64(), id, 0, state);
}

__device__ float PU_GPU(curandState *state) {
    //    Picks a random number from a uniform distribution

    float Random = curand_uniform(state);

    return Random;
}

__global__ void kernel(char *cells, double *vtemp, float *ut, float *ecl, float *inf, float *th, float *epl, float *ipl, float *rt, systemConstantsStruct constant, int cell2cell, int freecell, curandState *state, int NumberOfLayers, bool regensAllowed, float* timeDead) {
    int Row = threadIdx.x + blockIdx.x * blockDim.x;
    int Column = threadIdx.y + blockIdx.y * blockDim.y;

    int NX = (2 * NumberOfLayers - 1);
    int NY = (2 * NumberOfLayers - 1);
    int NXNY = NX * NY;

    if ((cells[Row + NX * Column + NXNY * 0] != 'o') && (Row + NX * Column + NXNY < 2 * NXNY)) {
        // Virus Spreads
        int AboveRow = Row - 1;        // row coordinate above cell
        int LeftColumn = Column - 1;   // column coordinate left of cell
        int BelowRow = Row + 1;        // row coordinate below cell
        int RightColumn = Column + 1;  // column coordinate right of cell

        //          if the cell one row up doesn't exist, it's taken out of the equation
        if (AboveRow < 0) {
            AboveRow = Row;
        }
        //          if the cell one column to the left doesn't exist, it's taken out of the equation
        if (LeftColumn < 0) {
            LeftColumn = Column;
        }
        //          if the cell one row down doesn't exist, it's taken out of the equation
        if (BelowRow > (NY - 1)) {
            BelowRow = Row;
        }
        //          if the cell one column to the right doesn't exist, it's taken out of the equation
        if (RightColumn > (NX - 1)) {
            RightColumn = Column;
        }

        if (cells[AboveRow + NX * Column + NXNY * 0] == 'o') {
            AboveRow = Row;
        }
        if (cells[AboveRow + NX * RightColumn + NXNY * 0] == 'o') {
            AboveRow = Row;
            RightColumn = Column;
        }
        if (cells[Row + NX * RightColumn + NXNY * 0] == 'o') {
            RightColumn = Column;
        }
        if (cells[BelowRow + NX * Column + NXNY * 0] == 'o') {
            BelowRow = Row;
        }
        if (cells[Row + NX * LeftColumn + NXNY * 0] == 'o') {
            LeftColumn = Column;
        }
        if (cells[BelowRow + NX * LeftColumn + NXNY * 0] == 'o') {
            BelowRow = Row;
            LeftColumn = Column;
        }

        double NNN = (vtemp[AboveRow + NX * Column + NXNY * 0] + vtemp[AboveRow + NX * RightColumn + NXNY * 0] + vtemp[Row + NX * RightColumn + NXNY * 0] + vtemp[BelowRow + NX * Column + NXNY * 0] + vtemp[Row + NX * LeftColumn + NXNY * 0] + vtemp[BelowRow + NX * LeftColumn + NXNY * 0]);

        double p = 0.0;
        if (cells[Row + NX * Column + NXNY * 0] == 'i') {
            p = constant.rho;
        }
        double VirusProduced = p * constant.timestep;
        double VirusDecay = constant.c * vtemp[Row + NX * Column + NXNY * 0] * constant.timestep;
        double VirusOut = 4.0 * constant.Dtsx2 * vtemp[Row + NX * Column + NXNY * 0];
        double VirusIn = 2.0 * constant.Dtsx2 * NNN / 3.0;

        __syncthreads();

        vtemp[Row + NX * Column + NXNY * 1] = vtemp[Row + NX * Column + NXNY * 0] + VirusProduced - VirusOut + VirusIn - VirusDecay;
        if (vtemp[Row + NX * Column + NXNY * 1] < pow(10.0, -10.0)) {
            vtemp[Row + NX * Column + NXNY * 1] = 0.0;
        }

        // The Cell behavior

        // Regen
        if (cells[Row + NX * Column + NXNY * 0] == 'd') {
            if (regensAllowed) {
                int index = Row + NX * Column;
                if (ut[index] > (inf[index] + ecl[index] + th[index] + rt[index])) {
                    th[index] = inf[index] + ecl[index] + th[index] + timeDead[index];
                    ecl[index] = epl[index];
                    inf[index] = ipl[index];
                    cells[Row + NX * Column + NX * NY * 1] = 'h';
                }
            }
        } 
        
        if (cells[Row + NX * Column + NXNY * 0] == 'i') {
            // Infectied
            if (ut[Row + NX * Column] > (inf[Row + NX * Column] + ecl[Row + NX * Column] + th[Row + NX * Column])) {
                cells[Row + NX * Column + NXNY * 1] = 'd';
                if (CODETESTINGCONDITIONS == 1) {
                    cells[Row + NX * Column + NXNY * 1] = 'i';
                }
            }
        } else if (cells[Row + NX * Column + NXNY * 0] == 'e') {
            // Eclipse
            if (ut[Row + NX * Column] > (ecl[Row + NX * Column] + th[Row + NX * Column])) {
                cells[Row + NX * Column + NXNY * 1] = 'i';
                inf[Row + NX * Column] = inf[Row + NX * Column] + ipl[Row + NX * Column];
            }
        } else if (cells[Row + NX * Column + NXNY * 0] == 'h') {
            // Healthy
            th[Row + NX * Column] = th[Row + NX * Column] + constant.timestep;

            if (cell2cell == 1) {
                // Cell to cell transmission
                int AboveRow = Row - 1;        // row coordinate above cell
                int LeftColumn = Column - 1;   // column coordinate left of cell
                int BelowRow = Row + 1;        // row coordinate below cell
                int RightColumn = Column + 1;  // column coordinate right of cell

                //        if the cell one row up doesn't exist, it's taken out of the equation
                if (AboveRow < 0) {
                    AboveRow = 0;
                }
                //        if the cell one column to the left doesn't exist, it's taken out of the equation
                if (LeftColumn < 0) {
                    LeftColumn = 0;
                }
                //        if the cell one row down doesn't exist, it's taken out of the equation
                if (BelowRow > NY - 1) {
                    BelowRow = 0;
                }
                //        if the cell one column to the right doesn't exist, it's taken out of the equation
                if (RightColumn > NX - 1) {
                    RightColumn = 0;
                }

                if (PU_GPU(state) < constant.probi * constant.timestep) {
                    if (cells[Row + NX * LeftColumn + NXNY * 0] == 'i') {
                        cells[Row + NX * Column + NXNY * 1] = 'e';
                    }
                    if (cells[Row + NX * RightColumn + NXNY * 0] == 'i') {
                        cells[Row + NX * Column + NXNY * 1] = 'e';
                    }
                    if (cells[AboveRow + NX * Column + NXNY * 0] == 'i') {
                        cells[Row + NX * Column + NXNY * 1] = 'e';
                    }
                    if (cells[BelowRow + NX * Column + NXNY * 0] == 'i') {
                        cells[Row + NX * Column + NXNY * 1] = 'e';
                    }
                    if (cells[AboveRow + NX * RightColumn + NXNY * 0] == 'i') {
                        cells[Row + NX * Column + NXNY * 1] = 'e';
                    }
                    if (cells[BelowRow + NX * LeftColumn + NXNY * 0] == 'i') {
                        cells[Row + NX * Column + NXNY * 1] = 'e';
                    }

                    ecl[Row + NX * Column] = epl[Row + NX * Column];
                }
            }

            if (freecell == 1) {
                // Cell free transmission
                double probablity = PU_GPU(state);
                double pnoinfect = exp(-vtemp[Row + NX * Column + NXNY * 0] * constant.beta * constant.timestep);
                if (probablity >= pnoinfect) {
                    cells[Row + NX * Column + NXNY * 1] = 'e';
                    ecl[Row + NX * Column] = epl[Row + NX * Column];
                }
            }
        }

        // The Universal Time for the cells is kept here (ut)
        ut[Row + NX * Column] = ut[Row + NX * Column] + constant.timestep;
        vtemp[Row + NX * Column + NXNY * 0] = vtemp[Row + NX * Column + NXNY * 1];
        cells[Row + NX * Column + NXNY * 0] = cells[Row + NX * Column + NXNY * 1];
    }
}

int main(int argc,char *argv[]) {
    // Checks for Heisenberg status of viral diffusion
    if (D * timestep / pow(deltxprime, 2.0) > 0.5) {
        printf("%.1f", D * timestep / pow(deltxprime, 2.0));
        printf("CHANGE PARAMETERS TO FIT DIFFUSION LIMITS. VALUE MUST BE UNDER 0.5. VALUE SHOWN ABOVE");
        exit(0);
    }
    // Clear Terminal
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

    // double MOI[15] = {pow(10, -1), 3 * pow(10, -1), 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25};
    double MOI[] = { pow(10, -2) };
    // double MOI[6] = {0.1, 1, 10, 15, 20, 25};
    // double MOI[21] = {5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25};
    // double MOI[5] = {70, 80, 100, 150, 250};
    // double MOI[5] = {100, 1000, 10000, 100000, 1000000};
    // double MOI[3] = {100, 1000, 10000};
    // double MOI[1] = {1000000};
    // float MOI[5] = {1,10,50,100,200};
    // float MOI[6] = {powf(10,0), powf(10,-1), powf(10,-2), powf(10,-3), powf(10,-4), powf(10,-5)};
    // float MOI[5] = {powf(10,-1), powf(10,-2), powf(10,-3), powf(10,-4), powf(10,-5)};
    // float MOI[3] = {powf(10,-3), powf(10,-4), powf(10,-5)};

    // float probi[4] = {0.2, 0.4, 0.6, 0.8};
    float probi[1] = {0.2};

    for (int q = 0; q < (sizeof(MOI) / sizeof(MOI[0])); q++) {
        for (int k = 0; k < (sizeof(probi) / sizeof(probi[0])); k++) {
            // Loop For The number Of Simulations To Run Per Setting
            for (int BigIndex = 0; BigIndex < NumberOfRuns; BigIndex++) {
                //        auto start = chrono::high_resolution_clock::now();

                //        printf("\nStarting run %d\n", (BigIndex+1));

                // Creating Save Path
                creatingPathToFolderAndDirectory(StartRuns + BigIndex, NumberOfLayers, MOI[q], regenParameter);
                // Creating placeholder variables for multipy runs
                int cell2cell = CELL2CELL;
                int freecell = FREECELL;

                // Building Cells
                creatingCellLocations();

                // Number of Cells
                // Number of initial infected cells
                int Ni = NumberOfCells * MOI[q];
                //        if(Ni < 1){ printf("Use larger MOI"); exit(0);}
                int Nx = (2 * NumberOfLayers - 1);  // Range of cells on x axis
                int Ny = (2 * NumberOfLayers - 1);  // Range of cells on y axis

                // Makeing empty matrices
                allocateMemory(Nx, Ny);

                // Initializing
                initailConditions(Nx, Ny, MOI[q]);

                //        //Deletes files and initial with values
                //        if(BigIndex == 0){
                //            printToFileCellAndVirusInitial(Nx, Ny, NumberOfLayers);
                //        }

                printToFileCellAndVirusAnalysisInitial(Nx, Ny);

                // Infects a random cell, now seen as (e)
                if (INITIALCELL == 1) {
                    infectANumberOfCellsRandomly(Nx, Ny, Ni);
                }

                if (RUNCPU == 0) {
                    cudaMalloc((void **)&state, Nx * Ny * sizeof(int));
                    errorCheck("cudaMalloc Random Setup");
                    cuRand_Setup<<<GridConfig, BlockConfig>>>(state);
                    errorCheck("Random Setup");

                    loadConstants(MOI[q], probi[k]);

                    deviceSetupAndMemoryAllocation(Nx, Ny);

                    cudaMemcpy(cells_GPU, cells, Nx * Ny * 2 * sizeof(char), cudaMemcpyHostToDevice);
                    errorCheck("cudaMemcpy cells HtoD");
                    cudaMemcpy(vtemp_GPU, vtemp, Nx * Ny * 2 * sizeof(double), cudaMemcpyHostToDevice);
                    errorCheck("cudaMemcpy vtemp HtoD");
                    cudaMemcpy(ut_GPU, ut, Nx * Ny * sizeof(float), cudaMemcpyHostToDevice);
                    errorCheck("cudaMemcpy ut HtoD");

                    cudaMemcpy(ecl_GPU, ecl, Nx * Ny * sizeof(float), cudaMemcpyHostToDevice);
                    errorCheck("cudaMemcpy ecl HtoD");
                    cudaMemcpy(inf_GPU, inf, Nx * Ny * sizeof(float), cudaMemcpyHostToDevice);
                    errorCheck("cudaMemcpy inf HtoD");
                    cudaMemcpy(th_GPU, th, Nx * Ny * sizeof(float), cudaMemcpyHostToDevice);
                    errorCheck("cudaMemcpy th HtoD");

                    cudaMemcpy(EclipsePhaseLength_GPU, EclipsePhaseLength, Nx * Ny * sizeof(float), cudaMemcpyHostToDevice);
                    errorCheck("cudaMemcpy EclipsePhaseLength HtoD");
                    cudaMemcpy(InfectionPhaseLength_GPU, InfectionPhaseLength, Nx * Ny * sizeof(float), cudaMemcpyHostToDevice);
                    errorCheck("cudaMemcpy InfectionPhaseLength HtoD");
                    
                    cudaMemcpy(RegenTime_GPU, RegenTime, Nx * Ny * sizeof(float), cudaMemcpyHostToDevice);
                    errorCheck("cudaMemcpy RegenTime HtoD");
                    cudaMemcpy(timeDead_GPU, timeDead, Nx * Ny * sizeof(float), cudaMemcpyHostToDevice);
                    errorCheck("cudaMemcpy timeDead HtoD");
                }

                // Runs simulation
                int NumberofTimeSteps = endtime / timestep;
                int NumberofSavedTimeSteps = NumberofTimeSteps / Save;
                float timestepcount = 0.0;  // equal to the number of ts elapsed
                while (timestepcount < (NumberofTimeSteps - 1)) {
                    if (RUNCPU == 0) {
                        kernel<<<GridConfig, BlockConfig>>>(cells_GPU, vtemp_GPU, ut_GPU, ecl_GPU, inf_GPU, th_GPU, EclipsePhaseLength_GPU, InfectionPhaseLength_GPU, RegenTime_GPU, SystemConstants, cell2cell, freecell, state, NumberOfLayers, regensAllowed, timeDead_GPU);
                    } else {
                        // Cerial Viral Transmission
                        //                cerialViralTransmission(Nx, Ny, cell2cell, freecell, probi[k]);
                        //                modifiedCerialViralTransmission(Nx, Ny, cell2cell, freecell, probi[k]);
                    }

                    if ((int(timestepcount) % Save) == 0) {
                        if (RUNCPU == 0) {
                            cudaMemcpy(cells, cells_GPU, Nx * Ny * 2 * sizeof(char), cudaMemcpyDeviceToHost);
                            errorCheck("cudaMemcpy cells DtoH");
                            cudaMemcpy(vtemp, vtemp_GPU, Nx * Ny * 2 * sizeof(double), cudaMemcpyDeviceToHost);
                            errorCheck("cudaMemcpy vtemp DtoH");
                        }

                        if (timestepcount != 0.0) {
                            printToFileCellAndVirusAnalysis(timestepcount * timestep);
                        }

                        // Dish analysis
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

                        //                //Prints status of cells  virus
                        //                if(BigIndex == 0){
                        //                    printToFileCellAndVirus(Nx, Ny, NumberOfLayers);
                        //                }
                    }

                    // Number of days completed
                    //            if((timestepcount%(24*int(1/timestep))) == 0){
                    //                printf("%.0f Day\n",(timestepcount*timestep)/24);
                    //            }

                    //            if((NumberInfected1 == 0) && (NumberEclipse1 == 0)){
                    //                cell2cell = 0;
                    //                freecell = 0;
                    //            }
                    //            else{
                    //                cell2cell = CELL2CELL;
                    //                freecell = FREECELL;
                    //            }
                    // End Code when all cells are not healthy
		if ((AmountOfVirus < pow(10, 1.0)) && (NumberDead1 == NumberOfCells)) {
                    break;
                }

                    timestepcount = timestepcount + 1.0;
                }

                // Writes a file with all of our parameters/variables
                createParameterFile(timestep, NumberofSavedTimeSteps, endtime, timestepcount, AmountOfVirus, rho, D, deltxprime, c, probi[k]);

                //        printf("\nMOI(%.1f) probi(%.1f): %d of %d Runs Done\n", log10(MOI[q]), probi[k], (StartRuns+BigIndex+1), NumberOfRuns);

                freeMemory();

                //        auto finish = std::chrono::high_resolution_clock::now();
                //        chrono::duration<double> elapsed = finish - start;
                //        cout << "Elapsed time: " << elapsed.count() << " s";
            }
        }
    }
    printf("\nPROGRAM DONE\n");
}
