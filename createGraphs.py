# Compatable with version 1.9 C Viral Trans

import numpy
import datetime
import matplotlib.pyplot as plt
import os
import sys
import datetime
import shutil

currentDT = datetime.datetime.now()

class Tc:
    # used for color printing
    INFO = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    END = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

print(Tc.BOLD, Tc.UNDERLINE, Tc.INFO, "Beginning Graph Analysis", Tc.END)

MOI = [-2.0] # can be multiple values of MOI
NumberOfRuns = 4
runFolders = 2
date1 = sys.argv[1]
date2 = "{}-{}".format(str(currentDT.month).zfill(2), str(currentDT.day).zfill(2)) # example: 03-14

# /run/media/asher/Backups/researchApprentice/scpTCU/  SCP-Mar20  /ahaun/Ccode/run1/ViralModel  /03-14/  193_4-FREECELL_-2.0-MOI_0.100-RP

OuterPath = r"/run/media/asher/Backups/researchApprentice/scpTCU/{}/ahaun/Ccode/ViralModel".format(date1)
NumberOfLayers = 193
regenParams = ["{:.3f}".format(0.003), "{:.3f}".format(0.010), "{:.3f}".format(0.030), "{:.3f}".format(0.100), "{:.3f}".format(0.400), "{:.3f}".format(1.000), "{:.3f}".format(3.000)] # must be floats
cellToCell = False

BigMOI = BigTimeArray = numpy.empty([len(MOI), NumberOfRuns],dtype=tuple )

BigTimeArray = numpy.empty([NumberOfRuns],dtype=tuple )
BigHealthyArray = numpy.empty([NumberOfRuns],dtype=tuple )
BigEclipseArray = numpy.empty([NumberOfRuns],dtype=tuple )
BigInfectedArray = numpy.empty([NumberOfRuns],dtype=tuple )
BigDeadArray = numpy.empty([NumberOfRuns],dtype=tuple )
BigVirusArray = numpy.empty([NumberOfRuns],dtype=tuple )

LongestTime = 0.0

# File structure looks like this (using example dates):
# SCP-Mar30/.../Ccode/
#   run1/
#       ViralModel/
#           03-04/
#               607_0-CELL2CELL_-2.0-MOI_0.1-RP/
#                   PerTimeStep.txt
#           03-05/
#               607_0-CELL2CELL_-2.0-MOI_0.1-RP/
#                   PerTimeStep.txt
#   run2/
#       ViralModel/
#           03-04/
#               607_0-FREECELL_-2.0-MOI_0.1-RP/
#                   PerTimeStep.txt
#           03-05/
#               607_0-FREECELL_-2.0-MOI_0.1-RP/
#                   PerTimeStep.txt

# We want it to look like this:
# SCP-Mar30/.../Ccode/
#   run1/
#       ViralModel/
#           607_0-CELL2CELL_-2.0-MOI_0.1-RP/
#               PerTimeStep.txt
#           607_0-CELL2CELL_-2.0-MOI_0.1-RP/
#               PerTimeStep.txt
#   run2/
#       ViralModel/
#           607_0-FREECELL_-2.0-MOI_0.1-RP/
#               PerTimeStep.txt
#           607_0-FREECELL_-2.0-MOI_0.1-RP/
#               PerTimeStep.txt

# this program also accounts for the seperate run#/ folders
# Note: folders that contain PerTimeStep.txt are going to be refered to as data folders (ie 607_0-CELL2CELL_-2.0-MOI_0.1-RP/)

for d in range(runFolders):
    # looks in ...../run1/ and ...../run2/
    for j in range(len(MOI)): # only using one MOI, this loop ONLY runs once.

        # calculating path to each viral model
        OUTERPATH_BEFORE = os.path.join(OuterPath).replace('\\','//')
        index = OUTERPATH_BEFORE.find('/ViralModel')

        runNumber = "/run{}".format(d+1)

        OUTERPATH = OUTERPATH_BEFORE[:index] + '{}'.format(runNumber) + OUTERPATH_BEFORE[index:]

        # Uses some BASH code to copy the data folders (example: 607_0-FREECELL_-2.0-MOI_0.1-RP/) up one directory into ViralModel/
        print(Tc.GREEN, "Path: ", OUTERPATH, "\n Run Folder: ", runNumber)
        print(Tc.BLUE, "Moving contents of subdirectories up one")
        os.system(
        """
        for d in {}/*; do
        cd "$d"
        cp -r * ../
        cd ..
        done
        ls
        """.format(OUTERPATH))

        print(Tc.END)

        # create Analysis folder for each run#/ folder
        directory = os.path.join("".join(OUTERPATH+ "_" +str(MOI[j])+ "_" +"Analysis"))
        if not os.path.exists(directory):
            os.makedirs(directory)

        # Finalizing the path to the run#/ViralModel/ for this index of 'd' in runFolders
        OUTER_Path_to_Folder = os.path.abspath(directory)

        # Actually beginning the graphing
        for k in regenParams:
            print("\n", Tc.INFO, "Regen Parameter: ", str(k))
            for i in range(NumberOfRuns):
                print(Tc.BLUE, " Run: ", i, Tc.END)

                # Figure out the name of the folder to look in for PerTimeStep.txt to extract data from
                if cellToCell:
                    InnerPath = r"/{}_{}-CELL2CELL_{}-MOI_{}-RP".format(NumberOfLayers, str(i), MOI[j], k)
                else:
                    InnerPath = r"/{}_{}-FREECELL_{}-MOI_{}-RP".format(NumberOfLayers, str(i), MOI[j], k)

                # Create a folder inside of InnerPath (each data folder) call "Analysis/", this will contain the graphs
                INNERPATH = os.path.join(OUTERPATH+InnerPath).replace('\\','//')
                directory = os.path.join("".join(INNERPATH + "/Analysis"))
                if not os.path.exists(directory):
                    os.makedirs(directory)
                Inner_Path_to_Folder = os.path.abspath(directory)

                # Begin harvesting the data from the data folder
                TimeArray = []
                HealthyArray = []
                EclipseArray = []
                InfectedArray = []
                DeadArray = []
                VirusArray = []

                with open(os.path.join(INNERPATH,"PerTimeStep.txt")) as CellInFile:
                    for cellline in CellInFile:
                        CellData = cellline
                        CellData = CellData.split(",")
                        del CellData[-1]

                        if(LongestTime < int(CellData[0])):
                            LongestTime = int(CellData[0])

                        TimeArray.append(int(CellData[0]))
                        HealthyArray.append(int(CellData[1]))
                        EclipseArray.append(int(CellData[2]))
                        InfectedArray.append(int(CellData[3]))
                        DeadArray.append(int(CellData[4]))
                        VirusArray.append(float(CellData[5]))

                # Used for calculating mean of data
                BigTimeArray[i] = TimeArray
                BigHealthyArray[i] = HealthyArray
                BigEclipseArray[i] = EclipseArray
                BigInfectedArray[i] = InfectedArray
                BigDeadArray[i] = DeadArray
                BigVirusArray[i] = VirusArray

                plt.rcParams['figure.figsize'] = [10, 10]
                fig, ax = plt.subplots()
                plt.plot(TimeArray, HealthyArray)
                plt.ylim(0, 10000) #
                plt.xticks(numpy.arange(0, TimeArray[-1], (TimeArray[-1]/10)))
                plt.close()
                fig.savefig(os.path.join(Inner_Path_to_Folder,"HealthyVsTime"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

                plt.rcParams['figure.figsize'] = [10, 10]
                fig, ax = plt.subplots()
                plt.plot(TimeArray, EclipseArray)
                plt.xticks(numpy.arange(0, TimeArray[-1], (TimeArray[-1]/10)))
                plt.close()
                fig.savefig(os.path.join(Inner_Path_to_Folder,"EclipseVsTime"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

                plt.rcParams['figure.figsize'] = [10, 10]
                fig, ax = plt.subplots()
                plt.plot(TimeArray, InfectedArray)
                plt.xticks(numpy.arange(0, TimeArray[-1], (TimeArray[-1]/10)))
                plt.close()
                fig.savefig(os.path.join(Inner_Path_to_Folder,"InfectedVsTime"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

                plt.rcParams['figure.figsize'] = [10, 10]
                fig, ax = plt.subplots()
                plt.plot(TimeArray, DeadArray)
                plt.xticks(numpy.arange(0, TimeArray[-1], (TimeArray[-1]/10)))
                plt.close()
                fig.savefig(os.path.join(Inner_Path_to_Folder,"DeadVsTime"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

                plt.rcParams['figure.figsize'] = [10, 10]
                fig, ax = plt.subplots()
                plt.yscale("log")
                plt.plot(TimeArray, VirusArray)
                plt.xticks(numpy.arange(0, TimeArray[-1], (TimeArray[-1]/10)))
                plt.title("Virus vs Time")
                plt.ylabel("Amount of Virus")
                plt.xlabel("Time-(hr)")
                plt.close()
                fig.savefig(os.path.join(Inner_Path_to_Folder,"VirusVsTime"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

                fo = open(os.path.join(Inner_Path_to_Folder,"data"+".txt"), "w")
                fo.write("Regen Parameter: {}".format(k))
                fo.write("Viral Peak: {}".format(numpy.max(VirusArray)))
                fo.write("Chronic Load: {}".format(numpy.average(VirusArray[150:400])))
                fo.write("Upslope: {}".format())
                fo.close()
                print("  Viral Peak: ", numpy.max(VirusArray))
                print("  Chronic Load: {}".format(numpy.average(VirusArray[500:600])))
                print("  Upslope: {}".format(numpy.log(VirusArray[50])/50))

                BigMOI[j][i] = max(BigVirusArray[i])

            # calculate the median data arrays
            MedianTimeArray = []
            MedianHealthyArray = []
            MedianEclipseArray = []
            MedianInfectedArray = []
            MedianDeadArray = []
            MedianVirusArray = []

            InbetweenTimeArray = []
            InbetwHealthyArray = []
            InbetwEclipseArray = []
            InbetwInfectedArray = []
            InbetwDeadArray = []
            InbetwVirusArray = []

            for i in range(LongestTime):
                for k in range(NumberOfRuns):
                    try:
                        InbetweenTimeArray.append(BigTimeArray[k][i])
                    except IndexError:
                        "Nothing"
                    try:
                        InbetwHealthyArray.append(BigHealthyArray[k][i])
                    except IndexError:
                        "Nothing"
                    try:
                        InbetwEclipseArray.append(BigEclipseArray[k][i])
                    except IndexError:
                        "Nothing"
                    try:
                        InbetwInfectedArray.append(BigInfectedArray[k][i])
                    except IndexError:
                        "Nothing"
                    try:
                        InbetwDeadArray.append(BigDeadArray[k][i])
                    except IndexError:
                        "Nothing"
                    try:
                        InbetwVirusArray.append(BigVirusArray[k][i])
                    except IndexError:
                        "Nothing"

                MedianTimeArray.append(numpy.median(InbetweenTimeArray))
                MedianHealthyArray.append(numpy.median(InbetwHealthyArray))
                MedianEclipseArray.append(numpy.median(InbetwEclipseArray))
                MedianInfectedArray.append(numpy.median(InbetwInfectedArray))
                MedianDeadArray.append(numpy.median(InbetwDeadArray))
                MedianVirusArray.append(numpy.median(InbetwVirusArray))

                InbetweenTimeArray.clear()
                InbetwHealthyArray.clear()
                InbetwEclipseArray.clear()
                InbetwInfectedArray.clear()
                InbetwDeadArray.clear()
                InbetwVirusArray.clear()

            plt.rcParams['figure.figsize'] = [10, 10]
            fig, ax = plt.subplots()
            plt.plot(MedianHealthyArray)
            #plt.ylim(0, 10000)

            plt.xticks(numpy.arange(0, LongestTime, LongestTime/10))
            plt.close()
            fig.savefig(os.path.join(OUTER_Path_to_Folder,"Median_HealthyVsTime{}".format(k) + ".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

            plt.rcParams['figure.figsize'] = [10, 10]
            fig, ax = plt.subplots()
            plt.plot(MedianEclipseArray)
            plt.xticks(numpy.arange(0, LongestTime, LongestTime/10))
            plt.close()
            fig.savefig(os.path.join(OUTER_Path_to_Folder,"Median_EclipseVsTime{}".format(k) + ".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

            plt.rcParams['figure.figsize'] = [10, 10]
            fig, ax = plt.subplots()
            plt.plot(MedianInfectedArray)
            plt.xticks(numpy.arange(0, LongestTime, LongestTime/10))
            plt.close()
            fig.savefig(os.path.join(OUTER_Path_to_Folder,"Median_InfectedVsTime{}".format(k) + ".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

            plt.rcParams['figure.figsize'] = [10, 10]
            fig, ax = plt.subplots()
            plt.plot(MedianDeadArray)
            plt.xticks(numpy.arange(0, LongestTime, LongestTime/10))
            plt.close()
            fig.savefig(os.path.join(OUTER_Path_to_Folder,"Median_DeadVsTime{}".format(k) + ".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

            plt.rcParams['figure.figsize'] = [10, 10]
            fig, ax = plt.subplots()
            plt.yscale("log")
            plt.plot(MedianVirusArray)
            plt.xticks(numpy.arange(0, LongestTime, LongestTime/10))
            plt.title("Median Virus vs Time")
            plt.ylabel("Median_ Virus Titer")
            plt.xlabel("Time-(hr)")
            plt.close()
            fig.savefig(os.path.join(OUTER_Path_to_Folder,"Median_VirusVsTime{}".format(k)+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

            plt.rcParams['figure.figsize'] = [10, 10]
            fig, ax = plt.subplots()
            for i in BigHealthyArray:
                plt.plot(i, color = "gray", alpha = 0.3)
            plt.plot(MedianHealthyArray, color = "red")
            plt.xticks(numpy.arange(0, LongestTime, LongestTime/10))
            plt.close()
            fig.savefig(os.path.join(OUTER_Path_to_Folder,"Median_HealthyVsTime-{}".format(k) + ".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

            plt.rcParams['figure.figsize'] = [10, 10]
            fig, ax = plt.subplots()
            for i in BigEclipseArray:
                plt.plot(i, color = "gray", alpha = 0.3)
            plt.plot(MedianEclipseArray, color = "red")
            plt.xticks(numpy.arange(0, LongestTime, LongestTime/10))
            plt.close()
            fig.savefig(os.path.join(OUTER_Path_to_Folder,"Median_EclipseVsTime-{}".format(k) + ".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

            plt.rcParams['figure.figsize'] = [10, 10]
            fig, ax = plt.subplots()
            for i in BigInfectedArray:
                plt.plot(i, color = "gray", alpha = 0.3)
            plt.plot(MedianInfectedArray, color = "red")
            plt.xticks(numpy.arange(0, LongestTime, LongestTime/10))
            plt.close()
            fig.savefig(os.path.join(OUTER_Path_to_Folder,"Median_InfectedVsTime-{}".format(k) + ".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

            plt.rcParams['figure.figsize'] = [10, 10]
            fig, ax = plt.subplots()
            for i in BigDeadArray:
                plt.plot(i, color = "gray", alpha = 0.3)
            plt.plot(MedianDeadArray, color = "red")
            plt.xticks(numpy.arange(0, LongestTime, LongestTime/10))
            plt.close()
            fig.savefig(os.path.join(OUTER_Path_to_Folder,"Median_DeadVsTime-{}".format(k) + ".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

            plt.rcParams['figure.figsize'] = [10, 10]
            fig, ax = plt.subplots()
            plt.yscale("log")
            for i in BigVirusArray:
                plt.plot(i, color = "gray", alpha = 0.3)
            plt.plot(MedianVirusArray, color = "red")
            plt.xticks(numpy.arange(0, LongestTime, LongestTime/10))
            plt.title("Median Virus vs Time")
            plt.ylabel("Median_ Virus Titer")
            plt.xlabel("Time-(hr)")
            plt.close()
            fig.savefig(os.path.join(OUTER_Path_to_Folder,"Median_VirusVsTime-{}".format(k) + ".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)


            plt.rcParams['figure.figsize'] = [10, 10]
            fig, ax = plt.subplots()
            histbins = numpy.linspace(min(BigMOI[j]),max(BigMOI[j]),20)
            hist, bins = numpy.histogram(BigMOI[j], bins=histbins)
            width = 0.7 * (bins[1] - bins[0])
            center = (bins[:-1] + bins[1:]) / 2
            ax.bar(center, hist, align='center', width=width)
            plt.close()
            fig.savefig(os.path.join(OUTER_Path_to_Folder,str(MOI[j])+"_Hist-{}".format(k) + ".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

print("Data Analysis Complete")
