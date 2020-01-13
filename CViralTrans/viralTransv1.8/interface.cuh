/**

  This interface must contain initializations for all variables and methods that are accessed globally/between files

*/

/** Mehods */
void deviceSetupAndMemoryAllocation(int Nx, int Ny);
float Te(float TauE, float ne);

float Ti(float TauI, float ni);

float PU1();

float exponentialDistro (double mean);

void errorCheck(const char *message);

void printInitialConditions(int SideLength, int RadiusOfCircle);

/** variables and constants **/
#define PI 3.1415926535897932f

#define CODETESTINGCONDITIONS 0
#define RUNCPU 1

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
extern float beta;
extern float rho;
extern float D;
extern float c;
extern float deltx;
extern float deltxprime;
extern float Dtsx2;

// Probability Constants
extern float TauI;
extern float TauE;
extern float ne;
extern float ni;
extern float probi;

extern double regenValues[];
extern int totalRegenerations;
extern double runParameter;

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
extern float* th_GPU;
extern float* ut;
extern float* ut_GPU;
extern float* EclipsePhaseLength;
extern float* EclipsePhaseLength_GPU;
extern float* InfectionPhaseLength;
extern float* InfectionPhaseLength_GPU;
extern int NumberOfCells;
extern int NumberDead;

extern int NumberDead1;
extern int NumberInfected1;
extern int NumberEclipse1;
extern int NumberHealthy1;
extern float AmountOfVirus;

extern curandState *state;
