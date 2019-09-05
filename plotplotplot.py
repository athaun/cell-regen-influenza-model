import numpy
import datetime
import matplotlib.pyplot as plt
import os

os.system('cls' if os.name=='nt' else 'clear')

MOI = [0.0, -1.0, -2.0, -3.0, -4.0, -5.0]
#MOI = [-5.0]
NumberOfRuns = 100
# C:\Users\Asher\Documents\Student_Research-master\Student_Research-master\ViralModel\7_0-FREECELL_0.0-MOI
OuterPath = r"/home/asher/Desktop/git/athaun/tcu/ViralModel"

BigMOI = BigTimeArray = numpy.empty([len(MOI), NumberOfRuns],dtype=tuple )

BigTimeArray = numpy.empty([NumberOfRuns],dtype=tuple )
BigHealthyArray = numpy.empty([NumberOfRuns],dtype=tuple )
BigEclipseArray = numpy.empty([NumberOfRuns],dtype=tuple )
BigInfectedArray = numpy.empty([NumberOfRuns],dtype=tuple )
BigDeadArray = numpy.empty([NumberOfRuns],dtype=tuple )
BigVirusArray = numpy.empty([NumberOfRuns],dtype=tuple )

LongestTime = 0.0

for j in range(len(MOI)):
    OUTERPATH = os.path.join(OuterPath).replace('\\','//')
    directory = os.path.join("".join(OUTERPATH+ "_" +str(MOI[j])+ "_" +"Analysis"))
    if not os.path.exists(directory):
        os.makedirs(directory)
    OUTER_Path_to_Folder = os.path.abspath(directory)

    for i in range(NumberOfRuns):

        os.system('cls' if os.name=='nt' else 'clear')
        print(str(MOI[j])+"\n")
        print(i)

        InnerPath = r"/7_"+str(i)+"-FREECELL_"+str(MOI[j])+"-MOI"

        INNERPATH = os.path.join(OuterPath+InnerPath).replace('\\','//')
        directory = os.path.join("".join(INNERPATH + "/Analysis"))
        if not os.path.exists(directory):
            os.makedirs(directory)
        Inner_Path_to_Folder = os.path.abspath(directory)

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

        BigTimeArray[i] = TimeArray
        BigHealthyArray[i] = HealthyArray
        BigEclipseArray[i] = EclipseArray
        BigInfectedArray[i] = InfectedArray
        BigDeadArray[i] = DeadArray
        BigVirusArray[i] = VirusArray

        plt.rcParams['figure.figsize'] = [10, 10]
        fig, ax = plt.subplots()
        plt.plot(TimeArray, HealthyArray)
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

        BigMOI[j][i] = max(BigVirusArray[i])

#        TimeArray.clear()
#        HealthyArray.clear()
#        EclipseArray.clear()
#        InfectedArray.clear()
#        DeadArray.clear()
#        VirusArray.clear()

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
    plt.xticks(numpy.arange(0, LongestTime, LongestTime/10))
    plt.close()
    fig.savefig(os.path.join(OUTER_Path_to_Folder,"Median_HealthyVsTime"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

    plt.rcParams['figure.figsize'] = [10, 10]
    fig, ax = plt.subplots()
    plt.plot(MedianEclipseArray)
    plt.xticks(numpy.arange(0, LongestTime, LongestTime/10))
    plt.close()
    fig.savefig(os.path.join(OUTER_Path_to_Folder,"Median_EclipseVsTime"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

    plt.rcParams['figure.figsize'] = [10, 10]
    fig, ax = plt.subplots()
    plt.plot(MedianInfectedArray)
    plt.xticks(numpy.arange(0, LongestTime, LongestTime/10))
    plt.close()
    fig.savefig(os.path.join(OUTER_Path_to_Folder,"Median_InfectedVsTime"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

    plt.rcParams['figure.figsize'] = [10, 10]
    fig, ax = plt.subplots()
    plt.plot(MedianDeadArray)
    plt.xticks(numpy.arange(0, LongestTime, LongestTime/10))
    plt.close()
    fig.savefig(os.path.join(OUTER_Path_to_Folder,"Median_DeadVsTime"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

    plt.rcParams['figure.figsize'] = [10, 10]
    fig, ax = plt.subplots()
    plt.yscale("log")
    plt.plot(MedianVirusArray)
    plt.xticks(numpy.arange(0, LongestTime, LongestTime/10))
    plt.title("Median Virus vs Time")
    plt.ylabel("Median_ Virus Titer")
    plt.xlabel("Time-(hr)")
    plt.close()
    fig.savefig(os.path.join(OUTER_Path_to_Folder,"Median_VirusVsTime"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

#    BigTimeArray[i] = TimeArray
#    BigHealthyArray[i] = HealthyArray
#    BigEclipseArray[i] = EclipseArray
#    BigInfectedArray[i] = InfectedArray
#    BigDeadArray[i] = DeadArray
#    BigVirusArray[i] = VirusArray

    plt.rcParams['figure.figsize'] = [10, 10]
    fig, ax = plt.subplots()
    for i in BigHealthyArray:
        plt.plot(i, color = "gray", alpha = 0.3)
    plt.plot(MedianHealthyArray, color = "red")
    plt.xticks(numpy.arange(0, LongestTime, LongestTime/10))
    plt.close()
    fig.savefig(os.path.join(OUTER_Path_to_Folder,"Median_HealthyVsTime"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

    plt.rcParams['figure.figsize'] = [10, 10]
    fig, ax = plt.subplots()
    for i in BigEclipseArray:
        plt.plot(i, color = "gray", alpha = 0.3)
    plt.plot(MedianEclipseArray, color = "red")
    plt.xticks(numpy.arange(0, LongestTime, LongestTime/10))
    plt.close()
    fig.savefig(os.path.join(OUTER_Path_to_Folder,"Median_EclipseVsTime"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

    plt.rcParams['figure.figsize'] = [10, 10]
    fig, ax = plt.subplots()
    for i in BigInfectedArray:
        plt.plot(i, color = "gray", alpha = 0.3)
    plt.plot(MedianInfectedArray, color = "red")
    plt.xticks(numpy.arange(0, LongestTime, LongestTime/10))
    plt.close()
    fig.savefig(os.path.join(OUTER_Path_to_Folder,"Median_InfectedVsTime"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

    plt.rcParams['figure.figsize'] = [10, 10]
    fig, ax = plt.subplots()
    for i in BigDeadArray:
        plt.plot(i, color = "gray", alpha = 0.3)
    plt.plot(MedianDeadArray, color = "red")
    plt.xticks(numpy.arange(0, LongestTime, LongestTime/10))
    plt.close()
    fig.savefig(os.path.join(OUTER_Path_to_Folder,"Median_DeadVsTime"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

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
    fig.savefig(os.path.join(OUTER_Path_to_Folder,"Median_VirusVsTime"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)


    plt.rcParams['figure.figsize'] = [10, 10]
    fig, ax = plt.subplots()
    histbins = numpy.linspace(min(BigMOI[j]),max(BigMOI[j]),20)
    hist, bins = numpy.histogram(BigMOI[j], bins=histbins)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    ax.bar(center, hist, align='center', width=width)
    plt.close()
    fig.savefig(os.path.join(OUTER_Path_to_Folder,str(MOI[j])+"_Hist"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

print("Data Analysis Complete")
