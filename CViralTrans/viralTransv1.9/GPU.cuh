#include "interface.cuh"

struct systemConstantsStruct {
    // NOTE  TO SELF: struct is basically the same as a dictionary in py, or class in java
    float MOI;
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

void loadConstants(float MOI) {
    SystemConstants.MOI = MOI;
    SystemConstants.beta = beta;
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

	GridConfig.x = (Nx - 1)/BlockConfig.x + 1;
	GridConfig.y = (Ny - 1)/BlockConfig.y + 1;
	GridConfig.z = 1;

	cudaMalloc((void**)&cells_GPU, Nx * Ny * 2 * sizeof(char));
	errorCheck("cudaMalloc cells Mem");
	cudaMalloc((void**)&vtemp_GPU, Nx * Ny * 2 * sizeof(float));
	errorCheck("cudaMalloc vtemp Mem");
	cudaMalloc((void**)&ut_GPU, Nx * Ny * sizeof(float));
	errorCheck("cudaMalloc ut Mem");

	cudaMalloc((void**)&ecl_GPU, Nx * Ny * sizeof(float));
	errorCheck("cudaMalloc ecl Mem");
	cudaMalloc((void**)&inf_GPU, Nx * Ny * sizeof(float));
	errorCheck("cudaMalloc inf Mem");
	cudaMalloc((void**)&th_GPU, Nx * Ny * sizeof(float));
	errorCheck("cudaMalloc th Mem");

	cudaMalloc((void**)&EclipsePhaseLength_GPU, Nx * Ny * sizeof(float));
	errorCheck("cudaMalloc EclipsePhaseLength Mem");
	cudaMalloc((void**)&InfectionPhaseLength_GPU, Nx * Ny * sizeof(float));
	errorCheck("cudaMalloc InfectionPhaseLength Mem");
}

__global__ void cuRand_Setup(curandState *state) {
    int Row = threadIdx.x + blockIdx.x * blockDim.x;
    int Column =  threadIdx.y + blockIdx.y * blockDim.y;
    int offsetx = blockDim.x * gridDim.x;

    int id = Row + offsetx * Column;
    curand_init (clock64(), id, 0, state);
}

__device__ float PU_GPU(curandState *state) {
    // Picks a random number from a uniform distribution
    float Random = curand_uniform(state);
    return Random;
}

__global__ void kernel(char *cells, float *vtemp, float *ut, float *ecl, float *inf, float *th,  float *epl, float *ipl, systemConstantsStruct constant, int cell2cell, int freecell, curandState *state, int NumberOfLayers) {

    int Row = threadIdx.x + blockIdx.x * blockDim.x;
    int Column =  threadIdx.y + blockIdx.y * blockDim.y;

    int NX = (2 * NumberOfLayers - 1);
    int NY = (2 * NumberOfLayers - 1);
    int NXNY = NX * NY;

    if ((cells[Row + NX * Column + NXNY * 0] != 'o') && (Row + NX * Column + NXNY < 2 * NXNY)) {
        // Virus Spreads
        int AboveRow = Row - 1;
        int LeftColumn = Column - 1;
        int BelowRow = Row + 1;
        int RightColumn = Column + 1;

        float rho2;
        if (cells[Row + NX * Column + NXNY * 0] == 'i') {
            rho2 = constant.rho;
        } else {
            rho2 = 0;
        }
        // where rho2 is a placeholder variable

        // if the cell one row up doesn't exist, it's taken out of the equation
        if (AboveRow < 0) {
            AboveRow = Row;
        }
        // if the cell one column to the left doesn't exist, it's taken out of the equation
        if (LeftColumn < 0) {
            LeftColumn = Column;
        }
        // if the cell one row down doesn't exist, it's taken out of the equation
        if (BelowRow > (NY - 1)) {
            BelowRow = Row;
        }
        // if the cell one column to the right doesn't exist, it's taken out of the equation
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

        float NNN = (vtemp[AboveRow + NX * Column + NXNY * 0] + vtemp[AboveRow + NX * RightColumn + NXNY * 0] + vtemp[Row + NX * RightColumn + NXNY * 0] + vtemp[BelowRow + NX * Column + NXNY * 0] + vtemp[Row + NX * LeftColumn + NXNY * 0] + vtemp[BelowRow + NX * LeftColumn + NXNY * 0]);

        float VirusProduced = rho2 * constant.timestep;
        float VirusDecay = constant.c * vtemp[Row + NX * Column + NXNY * 0]*constant.timestep;
        float VirusOut = 4.0 * constant.Dtsx2 * vtemp[Row + NX * Column + NXNY * 0];
        float VirusIn = 2.0 * constant.Dtsx2 * NNN / 3.0;

        __syncthreads();

        vtemp[Row + NX * Column + NXNY * 1] = vtemp[Row + NX * Column + NXNY * 0] + VirusProduced - VirusOut + VirusIn - VirusDecay;
        if (vtemp[Row + NX * Column + NXNY * 1] < pow(10.0,-10.0)) {
            vtemp[Row + NX * Column + NXNY * 1] = 0.0;
        }

        //The Cell behavior
        if (cells[Row + NX * Column + NXNY * 0] == 'i') {
            // Infectied
            if (ut[Row + NX * Column] > (inf[Row + NX * Column] + ecl[Row + NX * Column] + th[Row + NX * Column])) {
                cells[Row + NX * Column + NXNY * 1] = 'd';
                // if (CODETESTINGCONDITIONS == 1) {
                //     cells[Row + NX * Column + NXNY * 1] = 'i';
                // }
            }
        }
        else if (cells[Row + NX * Column + NXNY * 0] == 'e') {
            // Eclipse
            if (ut[Row + NX * Column] > (ecl[Row + NX * Column] + th[Row + NX * Column])) {
                cells[Row + NX * Column + NXNY * 1] = 'i';
                inf[Row + NX * Column] = inf[Row + NX * Column] + ipl[Row + NX * Column];
            }
        }
        else if (cells[Row + NX * Column + NXNY * 0] == 'h') {
            // Healthy
            th[Row + NX * Column] = th[Row + NX * Column] + constant.timestep;

            if (cell2cell == 1) {
                // Cell to cell transmission
                int AboveRow = Row - 1;   //row coordinate above cell
                int LeftColumn = Column - 1;   //column coordinate left of cell
                int BelowRow = Row + 1;   //row coordinate below cell
                int RightColumn = Column + 1;   //column coordinate right of cell

                // if the cell one row up doesn't exist, it's taken out of the equation
                if (AboveRow < 0) {
                    AboveRow = 0;
                }
                // if the cell one column to the left doesn't exist, it's taken out of the equation
                if (LeftColumn < 0) {
                    LeftColumn = 0;
                }
                // if the cell one row down doesn't exist, it's taken out of the equation
                if (BelowRow > NY - 1) {
                    BelowRow = 0;
                }
                // if the cell one column to the right doesn't exist, it's taken out of the equation
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
                float probablity = PU_GPU(state);
                float adaptedtimestep = constant.timestep; //variable time step
                float adaptedtimestepcount = 1.0;
                float pinfect = vtemp[Row + NX * Column + NXNY * 1] * constant.beta * adaptedtimestep;
                while(pinfect > 1.0) {
                    adaptedtimestep = adaptedtimestep / 2.0;
                    pinfect = vtemp[Row + NX * Column + NXNY * 1] * constant.beta * adaptedtimestep;
                    adaptedtimestepcount = adaptedtimestepcount * 2.0;
                }
                if (pinfect <= 1.0) {
                    if (adaptedtimestepcount != 1.0) {
                        pinfect = vtemp[Row + NX * Column + NXNY * 1] * constant.beta * adaptedtimestep;
                    }
                    while(adaptedtimestepcount != 1.0) {
                        if (probablity < pinfect) {
                            cells[Row + NX * Column + NXNY * 1] = 'e';
                            ecl[Row + NX * Column] = epl[Row + NX * Column];
                        }
                        adaptedtimestepcount = adaptedtimestepcount / 2.0;
                        adaptedtimestep = adaptedtimestep * 2.0;
                        pinfect = vtemp[Row + NX * Column + NXNY * 2] * constant.beta * adaptedtimestep;
                    }
                    if (adaptedtimestepcount == 1.0) {
                        vtemp[Row + NX * Column + NXNY * 1] = vtemp[Row + NX * Column + NXNY * 0] + VirusProduced - VirusOut + VirusIn - VirusDecay;
                        if (probablity < pinfect) {
                            cells[Row + NX * Column + NXNY * 1] = 'e';
                            ecl[Row + NX * Column] = epl[Row + NX * Column];
                        }
                    }
                }
            }
        }

        //The Universal Time for the cells is kept here (ut)
        ut[Row + NX * Column] = ut[Row + NX * Column] + constant.timestep;
        vtemp[Row + NX * Column + NXNY * 0] = vtemp[Row + NX * Column + NXNY * 1];
        cells[Row + NX * Column + NXNY * 0] = cells[Row + NX * Column + NXNY * 1];
    }
}
