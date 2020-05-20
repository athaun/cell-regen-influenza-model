#include "interface.cuh"

void runSimulation (float regenParameter) {
    for (int q = 0; q < (sizeof(MOI) / sizeof(MOI[0])); q ++) {
        for (int BigIndex = 0; BigIndex < NumberOfRuns; BigIndex ++) {
          // run the test as many times as necessary to get through all of the regen parameters, MOI values, and runs defined in NumberOfRuns

            printf("\nStarting run %d\n", (BigIndex + 1));


            //Creating Save Path
            creatingPathToFolderAndDirectory(StartRuns + BigIndex, NumberOfLayers, MOI[q], regenParameter);
            //Creating placeholder variables for multiple runs
            int cell2cell = CELL2CELL;
            int freecell = FREECELL;
            //Building Cells
            creatingCellLocations();

            std::cout << "Day\tHealthy\teclipse\tinfect\tdead\tregens\tViral Load" << std::endl;

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
                cudaMemcpy(RegenTime_GPU, RegenTime, Nx * Ny * sizeof(float), cudaMemcpyHostToDevice);
                errorCheck("cudaMemcpy RegenTime HtoD");
                cudaMemcpy(timeDead_GPU, timeDead, Nx * Ny * sizeof(float), cudaMemcpyHostToDevice);
                errorCheck("cudaMemcpy timeDead HtoD");
            }

            //Runs simulation
            int NumberofTimeSteps = endtime / timestep;
            int NumberofSavedTimeSteps = NumberofTimeSteps / Save;
            int timestepcount = 0;    //equal to the number of ts elapsed
            while(timestepcount < (NumberofTimeSteps - 1)) {

                if (RUNCPU == 0) {
                    // GPU code
                    kernel<<<GridConfig,BlockConfig>>>(cells_GPU, vtemp_GPU, ut_GPU, ecl_GPU, inf_GPU, th_GPU, EclipsePhaseLength_GPU, InfectionPhaseLength_GPU, RegenTime_GPU, SystemConstants, cell2cell, freecell, state, NumberOfLayers, regensAllowed, timeDead_GPU);
                } else {
                    CerialViralTransmission(Nx, Ny, cell2cell, freecell, regenParameter);
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

                    printToFileCellAndVirusAnalysis((timestepcount * timestep));
                }

                //Number of days completed
                if ((timestepcount % (24 * int(1 / timestep))) == 0) {
                    std::cout <<
                      (timestepcount * timestep) / 24 <<
                      "\t" << NumberHealthy <<
                      "\t" << NumberEclipse <<
                      "\t" << NumberInfected <<
                      "\t" << NumberDead <<
                      "\t" << totalRegenerations <<
                      "\t" << AmountOfVirus <<
                    std::endl;
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
