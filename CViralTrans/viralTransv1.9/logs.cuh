#include "interface.cuh"
/**
  this contains methods that print to log files, or manipulate the system
*/
void creatingPathToFolderAndDirectory (int BigIndex, int NumberOfLayers, float MOI, float regenParameter) {
    char TransmissionType[10] = "";
    if (CELL2CELL == 1) {
        if (FREECELL == 1) {
            strcat(TransmissionType,"Both");
        } else {
            strcat(TransmissionType,"CELL2CELL");
        }
    } else if (CELL2CELL == 0) {
	      if (FREECELL == 0) {
            strcat(TransmissionType,"Neither");
        } else {
            strcat(TransmissionType,"FREECELL");
        }
    }

    char Buffer[5]; // Buffer String For Conversion To Char
    char TheCurrentTime[50];
    time_t RawTime = time(NULL);
    tm* SpecificMoment = localtime(&RawTime);

    strcpy(Path_to_Folder, "");
    strcpy(directory, "");

    if (RUNCPU == 1) {
        strcat(Path_to_Folder,"ViralModel/");
    } else {
        strcat(Path_to_Folder,"/media / baylorfain / HDD / ViralModel/");
    }
    // strftime(TheCurrentTime, 50, "%m-%d/%I:%M", SpecificMoment);
    strftime(TheCurrentTime, 50, "%m-%d/", SpecificMoment);
    strcat(Path_to_Folder,TheCurrentTime);
    // strcat(Path_to_Folder,"_");
    sprintf(Buffer,"%d",NumberOfLayers);
    strcat(Path_to_Folder,Buffer);
    strcat(Path_to_Folder,"_");
    sprintf(Buffer,"%d",BigIndex);
    strcat(Path_to_Folder,Buffer);
    strcat(Path_to_Folder,"-");
    strcat(Path_to_Folder,TransmissionType);
    strcat(Path_to_Folder,"_");

    sprintf(Buffer,"%.1f",log10(MOI));
    strcat(Path_to_Folder,Buffer);
    strcat(Path_to_Folder,"-");
    strcat(Path_to_Folder,"MOI");

    strcat(Path_to_Folder, "_");
    sprintf(Buffer,"%.3f", regenParameter);
    strcat(Path_to_Folder, Buffer);
    strcat(Path_to_Folder, "-RP");
    // regenParameter

    strcat(directory,"mkdir -p ");
    strcat(directory,Path_to_Folder);
    int check = system(strdup(directory));
    if (check != 0) {
        exit(0);
    }
}

void printInitialConditions (int SideLength, int RadiusOfCircle, float regenParameter) {
  char File1[100] = "";
  strcat(File1,Path_to_Folder);
  strcat(File1,"/InitialCellLocations.txt");
  FILE *outfile1 = fopen(File1,"a");
  if (outfile1 == NULL) {
      printf("Error opening file!\n");
      exit(0);
  }

  for (int i = 0; i < ((2 * NumberOfLayers) - 1); i ++) {
      for (int j = 0; j < ((2 * NumberOfLayers) - 1); j ++) {
          fprintf(outfile1,"%c,",LocationData[i][j]);
      }
      fprintf(outfile1,"\n");
  }
  fclose(outfile1);

  char File2[100] = "";
  strcat(File2,Path_to_Folder);
  strcat(File2,"/Parameters.txt");
  FILE *outfile2 = fopen(File2, "w");
  if (outfile2 == NULL) {
      printf("Error opening file!\n");
      exit(0);
  }
  fprintf(outfile2, "Hexagon Side Length = %f\n", SideLength);
  fprintf(outfile2, "Number of Layers = %d\n", NumberOfLayers);
  fprintf(outfile2, "Radius of Circle = %f\n", RadiusOfCircle);
  fprintf(outfile2, "Number of Cells = %d\n", NumberOfCells);
  fclose(outfile2);

  printf("\033[95mHexagon Side Length = %f\n", SideLength);
  printf("Number of Layers = %d\n", NumberOfLayers);
  printf("Radius of Circle = %f\n", RadiusOfCircle);
  printf("Number of Cells = %d\n", NumberOfCells);
  printf("Exponential Distro = %f\n", regenParameter);
  printf("MOI = %d\n", MOI);
  printf("Saved to: %s\n\n\033[0m", Path_to_Folder);

}

void printToFileCellAndVirusInitial (int Nx, int Ny, int NumberOfLayers) {
    char File3[100] = "";
    strcat(File3, Path_to_Folder);
    strcat(File3, "/cells_over_time.txt");
    FILE *outfile3 = fopen(File3, "w");
    if (outfile3 == NULL) {
        printf("Error opening file!\n");
        exit(0);
    }
    for (int i = 0; i < ((2 * NumberOfLayers) - 1); i ++) {
        for (int j = 0; j < ((2 * NumberOfLayers) - 1); j ++) {
            fprintf(outfile3,"%c,",LocationData[i][j]);
        }
            fprintf(outfile3,"\n");
    }
    fclose(outfile3);

    char File4[100] = "";
    strcat(File4,Path_to_Folder);
    strcat(File4,"/virus_over_time.txt");
    FILE *outfile4 = fopen(File4,"w");
    if (outfile4 == NULL) {
        printf("Error opening file!\n");
        exit(0);
    }
    for (int i = 0; i < ((2 * NumberOfLayers) - 1); i ++) {
        for (int j = 0; j < ((2 * NumberOfLayers) - 1); j ++) {
            fprintf(outfile4,"%f,",0.0);
        }
            fprintf(outfile4,"\n");
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

void printToFileCellAndVirus (int Nx, int Ny, int NumberOfLayers) {
    char File5[100] = "";
    strcat(File5, Path_to_Folder);
    strcat(File5, "/cells_over_time.txt");
    FILE *outfile5 = fopen(File5, "a");
    if (outfile5 == NULL) {
        printf("Error opening file!\n");
        exit(0);
    }
    for (int i = 0; i<((2 * NumberOfLayers)-1); i++) {
        for (int j = 0; j<((2 * NumberOfLayers)-1); j++) {
            fprintf(outfile5,"%c,",cells[i + Nx * j + Nx * Ny * 0]);
        }
            fprintf(outfile5,"\n");
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
    for (int i = 0; i < ((2 * NumberOfLayers)-1); i++) {
        for (int j = 0; j < ((2 * NumberOfLayers)-1); j++) {
            fprintf(outfile6,"%f,",vtemp[i + Nx * j + Nx * Ny * 1]);
        }
        fprintf(outfile6,"\n");
    }
    fclose(outfile6);
}

void printToFileCellAndVirusAnalysis (float timestep) {
    char File8[100] = "";
    strcat(File8, Path_to_Folder);
    strcat(File8, "/PerTimeStep.txt");
    FILE *outfile8 = fopen(File8, "a");
    if (outfile8 == NULL) {
        printf("Error opening file!\n");
        exit(0);
    }
    /** why is it timestep + 1, and why wasn't it dvided by 24 */
    fprintf(outfile8,"%0.0f, %d, %d, %d, %d, %f, %d", timestep + 1, NumberHealthy1, NumberEclipse1, NumberInfected1, NumberDead1, AmountOfVirus, totalRegenerations);
    fprintf(outfile8,"\n");

    fclose(outfile8);
}

void createParameterFile (float timestep, int NumberofSavedTimeSteps, float endtime, float timestepcount, float AmountOfVirus, float rho, float D, float deltxprime, float c, float probi) {
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
    fprintf(outfile7, "D = %f\n", D);
    fprintf(outfile7, "delta x = %f\n", deltxprime);
    fprintf(outfile7, "c = %f\n", c);
    fprintf(outfile7, "Probability of cell to cell infection: %f\n", probi);
    fclose(outfile7);
}
