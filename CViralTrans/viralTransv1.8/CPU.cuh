#include "interface.cuh"
/**
  All code in this file is assumed to ONLY be used for CPU
*/

void creatingCellLocations() {
    float SideLength = (2.0 / 3.0);

    int RadiusScale = 0;
    int count = 0;
    for (int i = 0; i < NumberOfLayers; i ++) {
        RadiusScale ++;
        count ++;
    }
    float RadiusOfCircle = SideLength * RadiusScale;
    int NumberOfHexagons = count * 6 + 1;

    float** coord;
    int n = NumberOfHexagons;
    int m = 3;
    coord = (float**) calloc(n, sizeof(float*));
    for (int i = 0; i < n; i ++) {
       coord[i] = (float*) calloc(m, sizeof(float));
    }

    float** percyclecoord;
    n = NumberOfHexagons;
    m = 3;
    percyclecoord = (float**) calloc(n, sizeof(float*));
    for (int i = 0; i < n; i++) {
       percyclecoord[i] = (float*) calloc(m, sizeof(float));
    }

    int temp;
    for (int j = 0; j < NumberOfLayers; j ++) {
        for (int i = 0; i < (2 * j); i ++) {
            if (i < j) {
                temp = i;
            }
            percyclecoord[i + (j - 1) * j + 1][0] =  -temp - 1;
            percyclecoord[i + (j - 1) * j + 1][1] =   temp + j - i;
            percyclecoord[i + (j - 1) * j + 1][2] =  -j + 1 + i;
        }
    }
    float c0[3] = {percyclecoord[0][0], percyclecoord[0][1], percyclecoord[0][2]};
    coord[0][0] = c0[0];
    coord[0][1] = c0[1];
    coord[0][2] = c0[2];

    count = 0;
    for (int j = 0; j < (NumberOfHexagons / 3); j ++) {
        for (int i = 0; i<3; i++) {
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

    for (int i = 0; i < NumberOfHexagons; i ++) {
        xcoord = coord[i][0];
        if (coord[i][0] < xmin) {
            xmin = coord[i][0];
        }
        ycoord = (2.0 * sin(PI*(60.0 / 180.0)) * (coord[i][1] - coord[i][2]) / 3.0) + vi;
        dist = sqrtf(pow(double(xcoord - hi), 2.0) + pow(double(ycoord - vi), 2.0));
        if (dist >= RadiusOfCircle) {
            coord[i][0] = 5000.0;
            coord[i][1] = 0.0;
            coord[i][2] = 0.0;
        }
    }

    n = (2 * NumberOfLayers) - 1;
    m = (2 * NumberOfLayers) - 1;
    LocationData = (char**) malloc(n * sizeof(char*));
    for (int j = 0; j < n; j ++) {
        LocationData[j] = (char*) malloc(m * sizeof(char));
        for (int i = 0; i < m; i ++) {
            LocationData[j][i] = 'o';
       }
    }

    NumberOfCells = 0;
    for (int i = 0; i < NumberOfHexagons; i ++) {
        if (coord[i][0] != 5000.0) {
            LocationData[int(coord[i][2]) - int(xmin)][int(coord[i][0]) - int(xmin)] = 'h';
            NumberOfCells ++;
        }
    }

    printInitialConditions(SideLength, RadiusOfCircle);

    for (int i = 0; i < NumberOfHexagons; i++) {
       free(coord[i]);
       free(percyclecoord[i]);
    }
    free(coord);
    free(percyclecoord);
}

void allocateMemory(int Nx, int Ny) {
    cells = (char*) malloc(Nx * Ny * 2 * sizeof(char)); // Produces a matrix for the cells
    vtemp = (float*) malloc(Nx * Ny * 2 * sizeof(float)); // Produces a matrix that will track the amount virus above each cell
    ut = (float*) malloc(Nx * Ny * sizeof(float)); // Produces a univeral time matrix (ut)
    ecl = (float*) malloc(Nx * Ny * sizeof(float)); // Produces a time matrix for after eclipse phase (e)
    inf = (float*) malloc(Nx * Ny * sizeof(float)); // Produces a time matrix for after infection phase (i)
    th = (float*) malloc(Nx * Ny * sizeof(float)); // Produces a time matrix hor healthy cells (t)
    timeDead = (float*) malloc(Nx * Ny * sizeof(float));
    EclipsePhaseLength = (float*) malloc(Nx * Ny * sizeof(float)); // Produces an array of eclipse phase durations for cells
    InfectionPhaseLength = (float*) malloc(Nx * Ny * sizeof(float)); // Produces an array of infection phase durations for cells

}

void initailConditions(int Nx, int Ny) {
    for (int j = 0; j < Ny; j ++) {
        for (int i = 0; i < Nx; i ++) {
            for (int k = 0; k < 2; k ++) {
                cells[i + Nx * j + Nx * Ny * k] = LocationData[i][j];
                vtemp[i + Nx * j + Nx * Ny * k] = 0.0;
            }
            ut[i + Nx * j] = 0.0;
            ecl[i + Nx * j] = 0.0;
            inf[i + Nx * j] = 0.0;
            th[i + Nx * j] = 0.0;
            timeDead[i + Nx * j] = 0.0;
            EclipsePhaseLength[i + Nx * j] = Te(TauE,ne);
            InfectionPhaseLength[i + Nx * j]  = Ti(TauI,ni);
       }
    }
}

void infectRandomCells(int Nx, int Ny, int Ni) {
    if (CODETESTINGCONDITIONS == 1) {
        cells[(NumberOfLayers - 1) + Nx * (NumberOfLayers - 1) + Nx * Ny * 0] = 'i';
        cells[(NumberOfLayers - 1) + Nx * (NumberOfLayers - 1) + Nx * Ny * 1] = 'i'; //Only the center cell
    }
    else {
        srand(time(NULL));
        int randx, randy;
        int NumberOfInfectedCellsCount = 0;
        while(NumberOfInfectedCellsCount < Ni) {
            randx = (rand()%Nx);
            randy = (rand()%Ny);
            if ((cells[randx + Nx * randy + Nx * Ny * 0] != 'o') && (cells[randx + Nx * randy + Nx * Ny * 0] == 'h')) {
                cells[randx + Nx * randy + Nx * Ny * 0] = 'e';
                cells[randx + Nx * randy + Nx * Ny * 1] = 'e';
                NumberOfInfectedCellsCount ++;
            }
        }
    }
}

void modifiedCerialViralTransmission(int Nx, int Ny, int cell2cell, int freecell) {

        float te = Te(TauE,ne);
        float ti = Ti(TauI, ni);

        int NumberHealthy = 0;
        int NumberEclipse = 0;
        int NumberInfected = 0;
        int NumberVirus = 0;

        for (int y = 0; y < Ny; y ++) {
            for (int x = 0; x < Nx; x ++) {
                switch (cells[x + Nx * y + Nx * Ny * 0]) { // why multiply by 0?
                  case 'h':
                    NumberHealthy ++;
                    break;
                  case 'e':
                    NumberEclipse ++;
                    break;
                  case 'i':
                    NumberInfected ++;
                    break;
                  case 'o':
                    NumberVirus ++;
                    break;
                }
            }
        }

        int** LocationHealthy = (int**) malloc(NumberHealthy * sizeof(int*));
        int** LocationEclipse = (int**) malloc(NumberEclipse * sizeof(int*));
        int** LocationInfected = (int**) malloc(NumberInfected * sizeof(int*));
        int** LocationVirus = (int**) malloc(NumberVirus * sizeof(int*));

        for (int i = 0; i < NumberHealthy; i ++) {
           LocationHealthy[i] = (int*) malloc(2 * sizeof(int));
        }
        for (int i = 0; i < NumberEclipse; i ++) {
           LocationEclipse[i] = (int*) malloc(2 * sizeof(int));
        }
        for (int i = 0; i < NumberInfected; i ++) {
           LocationInfected[i] = (int*) malloc(2 * sizeof(int));
        }
        for (int i = 0; i < NumberVirus; i ++) {
           LocationVirus[i] = (int*) malloc(2 * sizeof(int));
        }

        int IndexerH = 0, IndexerE = 0, IndexerI = 0, IndexerO = 0;
        for (int y = 0; y < Ny; y ++) {
            for (int x = 0; x < Nx; x ++) {
                switch (cells[x + Nx * y + Nx * Ny * 0]) {
                  case 'h':
                    LocationHealthy[IndexerH][0] = x;
                    LocationHealthy[IndexerH][1] = y;
                    IndexerH ++;
                    break;
                  case 'e':
                    LocationEclipse[IndexerE][0] = x;
                    LocationEclipse[IndexerE][1] = y;
                    IndexerE ++;
                    break;
                  case 'i':
                    LocationInfected[IndexerI][0] = x;
                    LocationInfected[IndexerI][1] = y;
                    IndexerI ++;
                    break;
                  case 'o':
                    LocationVirus[IndexerO][0] = x;
                    LocationVirus[IndexerO][1] = y;
                    IndexerO ++;
                    break;
                }
            }
        }

        /**
          HealthyCells time
        */

        if (NumberHealthy != 0) {
            int Row, Column;
            for (int j = 0; j < NumberHealthy; j ++) {
                Row = LocationHealthy[j][0];
                Column = LocationHealthy[j][1];
                th[Row + Nx * Column] += timestep;
                // "th" is the time matrix for healthy cells
            }
        }

        /**
          Eclipse phase -> Infection
        */

        if (NumberEclipse != 0) {
            int Row, Column;
            for (int j = 0; j < NumberEclipse; j ++) {
                Row = LocationEclipse[j][0];
                Column = LocationEclipse[j][1];

                if ((ecl[Row + Nx * Column] + th[Row + Nx * Column]) < ut[Row + Nx * Column]) {
                    cells[Row + Nx * Column + Nx * Ny * 1] = 'i';
                    inf[Row + Nx * Column] = inf[Row + Nx * Column] + ti; // Ti(TauI, ni);
                }
            }
        }

        /**
          Infection spreads ================================================================================================
        */

        if (cell2cell == 1) {
            if (NumberInfected != 0) {
                int Row, Column;
                for (int j = 0; j < NumberInfected; j ++) {
                    Row = LocationInfected[j][0];
                    Column = LocationInfected[j][1];

                    // Binary toggles
                    int AboveRowExists = 1;
                    int LeftColumnExists = 1;
                    int BelowRowExists = 1;
                    int RightColumnExists = 1;

                    int AboveRow = Row - 1;        // row coordinate above cell
                    int LeftColumn = Column - 1;   // column coordinate left of cell
                    int BelowRow = Row + 1;        // row coordinate below cell
                    int RightColumn = Column + 1;  // column coordinate right of cell

                    // if the cell one row up doesn't exist, it's taken out of the equation
                    if (AboveRow < 0) {
                        AboveRowExists = 0;
                        AboveRow = 0;
                    }
                    // if the cell one column to the left doesn't exist, it's taken out of the equation
                    if (LeftColumn < 0) {
                        LeftColumnExists = 0;
                        LeftColumn = 0;
                    }
                    // if the cell one row down doesn't exist, it's taken out of the equation
                    if (BelowRow > Ny - 1) {
                        BelowRowExists = 0;
                        BelowRow = 0;
                    }
                    // if the cell one column to the right doesn't exist, it's taken out of the equation
                    if (RightColumn > Nx - 1) {
                        RightColumnExists = 0;
                        RightColumn = 0;
                    }

                    if (PU1() < probi * timestep) {
                        if ((LeftColumnExists == 1) && (cells[Row + Nx * LeftColumn + Nx * Ny * 0] != 'o')) {
                            if (cells[Row + Nx * LeftColumn + Nx * Ny * 0] == 'h') {
                                cells[Row + Nx * LeftColumn + Nx * Ny * 1] = 'e';
                                ecl[Row + Nx * LeftColumn] = te;
                            }
                        }

                        if ((RightColumnExists == 1) && (cells[Row + Nx * RightColumn + Nx * Ny * 0] != 'o')) {
                            if (cells[Row + Nx * RightColumn + Nx * Ny * 0] == 'h') {
                                cells[Row + Nx * RightColumn + Nx * Ny * 1] = 'e';
                                ecl[Row + Nx * RightColumn] = te;
                            }
                        }

                        if ((AboveRowExists == 1) && (cells[AboveRow + Nx * Column + Nx * Ny * 0] != 'o')) {
                            if (cells[AboveRow + Nx * Column + Nx * Ny * 0] == 'h') {
                                cells[AboveRow + Nx * Column + Nx * Ny * 1] = 'e';
                                ecl[AboveRow + Nx * Column] = te;
                            }
                        }

                        if ((BelowRowExists == 1) && (cells[BelowRow + Nx * Column + Nx * Ny * 0] != 'o')) {
                            if (cells[BelowRow + Nx * Column + Nx * Ny * 0] == 'h') {
                                cells[BelowRow + Nx * Column + Nx * Ny * 1] = 'e';
                                ecl[BelowRow + Nx * Column] = te;
                            }
                        }

                        if ((AboveRowExists == 1) && (RightColumnExists == 1) && (cells[AboveRow + Nx * RightColumn + Nx * Ny * 0] != 'o')) {
                            if (cells[AboveRow + Nx * RightColumn + Nx * Ny * 0] == 'h') {
                                cells[AboveRow + Nx * RightColumn + Nx * Ny * 1] = 'e';
                                ecl[AboveRow + Nx * RightColumn] = te;
                            }
                        }

                        if ((BelowRowExists == 1) && (LeftColumnExists == 1) && (cells[BelowRow + Nx * LeftColumn + Nx * Ny * 0] != 'o')) {
                            if (cells[BelowRow + Nx * LeftColumn + Nx * Ny * 0] == 'h') {
                                cells[BelowRow + Nx * LeftColumn + Nx * Ny * 1] = 'e';
                                ecl[BelowRow + Nx * LeftColumn] = te;
                            }
                        }
                    }
                }
            }
        }

        /**
          Virus Spreads
        */

        int Row, Column;
        for (int j = 0; j < NumberVirus; j ++) {
            Row = LocationVirus[j][0];
            Column = LocationVirus[j][1];

            int AboveRow = Row - 1;
            int LeftColumn = Column - 1;
            int BelowRow = Row + 1;
            int RightColumn = Column + 1;

            float rho2;
            if (cells[Row + Nx * Column + Nx * Ny * 0] == 'i') {
                rho2 = rho;
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
            if (BelowRow > (Ny - 1)) {
                BelowRow = Row;
            }
            // if the cell one column to the right doesn't exist, it's taken out of the equation
            if (RightColumn > (Nx - 1)) {
                RightColumn = Column;
            }

            if (cells[AboveRow + Nx * Column + Nx * Ny * 0] == 'o') {
                AboveRow = Row;
            }
            if (cells[AboveRow + Nx * RightColumn + Nx * Ny * 0] == 'o') {
                AboveRow = Row;
                RightColumn = Column;
            }
            if (cells[Row + Nx * RightColumn + Nx * Ny * 0] == 'o') {
                RightColumn = Column;
            }
            if (cells[BelowRow + Nx * Column + Nx * Ny * 0] == 'o') {
                BelowRow = Row;
            }
            if (cells[Row + Nx * LeftColumn + Nx * Ny * 0] == 'o') {
                LeftColumn = Column;
            }
            if (cells[BelowRow + Nx * LeftColumn + Nx * Ny * 0] == 'o') {
                BelowRow = Row;
                LeftColumn = Column;
            }

            float NNN = (vtemp[AboveRow + Nx * Column + Nx * Ny * 0] + vtemp[AboveRow + Nx * RightColumn + Nx * Ny * 0] + vtemp[Row + Nx * RightColumn + Nx * Ny * 0] + vtemp[BelowRow + Nx * Column + Nx * Ny * 0] + vtemp[Row + Nx * LeftColumn + Nx * Ny * 0] + vtemp[BelowRow + Nx * LeftColumn + Nx * Ny * 0]);

            float VirusProduced = rho2 * timestep;
            float VirusDecay = c * vtemp[Row + Nx * Column + Nx * Ny * 0]*timestep;
            float VirusOut = 4.0 * Dtsx2 * vtemp[Row + Nx * Column + Nx * Ny * 0];
            float VirusIn = 2.0 * Dtsx2 * NNN / 3.0;

            vtemp[Row + Nx * Column + Nx * Ny * 1] = vtemp[Row + Nx * Column + Nx * Ny * 0] + VirusProduced - VirusOut + VirusIn - VirusDecay;
            if (vtemp[Row + Nx * Column + Nx * Ny * 1] < pow(10.0,-10.0)) {
                vtemp[Row + Nx * Column + Nx * Ny * 1] = 0.0;
            }
            // Probability of infect adaptive time step
            if (freecell == 1) {
                float probability = PU1();
                float adaptedtimestep = timestep; // Variable time step
                float adaptedtimestepcount = 1.0;
                float pinfect = vtemp[Row + Nx * Column + Nx * Ny * 1] * beta * adaptedtimestep;
                while(pinfect > 1.0) {
                    adaptedtimestep /= 2.0;
                    pinfect = vtemp[Row + Nx * Column + Nx * Ny * 1] * beta * adaptedtimestep;
                    adaptedtimestepcount *= 2.0;
                }
                if (pinfect <= 1.0) {
                    if (adaptedtimestepcount != 1.0) {
                        pinfect = vtemp[Row + Nx * Column + Nx * Ny * 1]*beta * adaptedtimestep;
                    }
                    while(adaptedtimestepcount != 1.0) {
                        if (probability < pinfect) {
                            if (cells[Row + Nx * Column + Nx * Ny * 0] == 'h') {
                                cells[Row + Nx * Column + Nx * Ny * 1] = 'e';
                                ecl[Row + Nx * Column] = te;
                            }
                        }
                        adaptedtimestepcount /= 2.0;
                        adaptedtimestep *= 2.0;
                        pinfect = vtemp[Row + Nx * Column + Nx * Ny * 1]*beta * adaptedtimestep;
                    }
                    if (adaptedtimestepcount == 1.0) {
                        vtemp[Row + Nx * Column + Nx * Ny * 1] = vtemp[Row + Nx * Column + Nx * Ny * 0] + VirusProduced - VirusOut + VirusIn - VirusDecay;
                        if (probability < pinfect) {
                            if (cells[Row + Nx * Column + Nx * Ny * 0] == 'h') {
                                cells[Row + Nx * Column + Nx * Ny * 1] = 'e';
                                ecl[Row + Nx * Column] = te;
                            }
                        }
                    }
                }
            }
        }

        /**
          kills cells
        */
        if (NumberInfected != 0) {
            int Row;
            int Column;
            for (int j = 0; j < NumberInfected; j++) {
                Row = LocationInfected[j][0];
                Column = LocationInfected[j][1];
                if (ut[Row + Nx * Column] > (inf[Row + Nx * Column] + ecl[Row + Nx * Column] + th[Row + Nx * Column])) {
                    cells[Row + Nx * Column + Nx * Ny * 1] = 'd';
                    if (CODETESTINGCONDITIONS == 1) { // What is this used for?
                        cells[Row + Nx * Column + Nx * Ny * 1] = 'i';
                    }
                }
            }
        }

        for (int i = 0; i < NumberHealthy; i++) {
           free(LocationHealthy[i]);
        }
        for (int i = 0; i < NumberEclipse; i++) {
           free(LocationEclipse[i]);
        }
        for (int i = 0; i < NumberInfected; i++) {
           free(LocationInfected[i]);
        }
        for (int i = 0; i < NumberVirus; i++) {
           free(LocationVirus[i]);
        }
        free(LocationHealthy);
        free(LocationEclipse);
        free(LocationInfected);
        free(LocationVirus);

        for (int j = 0; j<Ny; j++) {
            for (int i = 0; i<Nx; i++) {
                vtemp[i + Nx * j + Nx * Ny * 0] = vtemp[i + Nx * j + Nx * Ny * 1];
                cells[i + Nx * j + Nx * Ny * 0] = cells[i + Nx * j + Nx * Ny * 1];
            }
        }

        //The Universal Time for the cells is kept here (ut)
        for (int j = 0; j < Ny; j ++) {
            for (int i = 0; i < Nx; i ++) {
                ut[i + Nx * j] += timestep;
            }
        }
}
