import numpy
import datetime
import matplotlib.pyplot as plt
import os
from scipy.stats import norm
from scipy.stats import shapiro
from scipy.stats import anderson
from scipy.stats import normaltest
from scipy.stats import ttest_1samp
import scipy.integrate as integrate
import matplotlib.ticker as ticker

# os.system('cls' if os.name=='nt' else 'clear')

def growthfit(TimeArray, VirusArray, IndexOfPeakTime):
    import scipy.optimize as optim
    import numpy as np
    if IndexOfPeakTime <= 2:
        b = 0.0
        perr = [0.0, 0.0, 0.0]
        r2 = 0.0
    else:
        def logistic(t, a, b, c):
            return a/(1+np.exp(-b*(t-c)))
        beginningtime = IndexOfPeakTime*5//10
        x = TimeArray[beginningtime:IndexOfPeakTime]
        y = np.divide(VirusArray[beginningtime:IndexOfPeakTime],VirusArray[IndexOfPeakTime])
        bounds = (0, [2,1,IndexOfPeakTime])
        (a,b,c),cov = optim.curve_fit(logistic, x, y, bounds=bounds)
        perr = np.sqrt(np.diag(cov))
        y_fit = numpy.empty([len(x)],dtype=float )
        for m in range(len(x)):
            y_fit[m] = logistic(x[m], a, b, c)
        # residual sum of squares
        ss_res = np.sum((y - y_fit) ** 2)
        # total sum of squares
        ss_tot = np.sum((y - np.mean(y)) ** 2)
        # r-squared
        r2 = 1 - (ss_res / ss_tot)
    return (b, perr[1], r2)

def decayfit(TimeArray,EclipseArray,InfectedArray,DeadArray,VirusArray):
    import scipy.optimize as optim
    import numpy as np
    notif = 0
    for ii in range(len(TimeArray)):
        if (EclipseArray[ii]==0)&(InfectedArray[ii]==0)&(DeadArray[ii]!=0)&(notif==0):
            StartTime = TimeArray[ii]
            notif = 1
        if (VirusArray[ii] <= 10**2)&(TimeArray[ii] > 5):
            EndTime = TimeArray[ii]
            break

    x = TimeArray[StartTime:EndTime]
    y = np.log(VirusArray[StartTime:EndTime])

    def linear(x, m, b):
        return m*x+b

    (m,b),cov = optim.curve_fit(linear, x, y)
    perr = np.sqrt(np.diag(cov))

    y_fit = numpy.empty([len(x)],dtype=float )
    for iii in range(len(x)):
        y_fit[iii] = linear(iii,m,b)
    # residual sum of squares
    ss_res = np.sum((y - y_fit) ** 2)
    # total sum of squares
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    # r-squared
    r2 = 1 - (ss_res / ss_tot)
    return(np.abs(m), perr[1], r2)

params = [
        0.003,
        0.01,
        0.03,
        0.1,
        0.3,
        0.6,
        0.8,
        1.0,
        1.3,
        1.7,
        2.0,
        3.0
        ]
#regenParams = ["{:.3f}".format(0.010), "{:.3f}".format(0.030), "{:.3f}".format(0.100), "{:.3f}".format(0.200), "{:.3f}".format(0.300), "{:.3f}".format(0.500),"{:.3f}".format(0.700), "{:.3f}".format(1.000), "{:.3f}".format(2.000), "{:.3f}".format(3.000)]
regenParams = []

for i in params:
    regenParams.append("{:.3f}".format(i))

# regenParams = ["{:.3f}".format(1.000)]
regenParamsFloat = []
for i in range(len(regenParams)):
    regenParamsFloat.append(float(regenParams[i]))
NumberOfRuns = 10

OuterPath = r"/media/asher/Backups/finalRun-7-14-20/free/ViralModel"
#OuterPath = r"/media/baylor/My Passport/BaylorFain/SARS-CoV-2_Runs/SARS-CoV-2"
#OuterPath = r"/home/baylor/Documents/Research/Coding/Viral_Transmission/ViralModel/02-03"


BigPeakVirus = numpy.empty([len(regenParams), NumberOfRuns],dtype=float )
BigTimeOfPeak = numpy.empty([len(regenParams), NumberOfRuns],dtype=float )
BigUpSlope = numpy.empty([len(regenParams), NumberOfRuns],dtype=float )
BigDownSlope = numpy.empty([len(regenParams), NumberOfRuns],dtype=float )
BigEndTime = numpy.empty([len(regenParams), NumberOfRuns],dtype=float )
BigAUC = numpy.empty([len(regenParams), NumberOfRuns],dtype=float )
BigChronicLoad = numpy.empty([len(regenParams), NumberOfRuns],dtype=float )

BigTimeArray = numpy.empty([NumberOfRuns],dtype=tuple )
BigHealthyArray = numpy.empty([NumberOfRuns],dtype=tuple )
BigEclipseArray = numpy.empty([NumberOfRuns],dtype=tuple )
BigInfectedArray = numpy.empty([NumberOfRuns],dtype=tuple )
BigDeadArray = numpy.empty([NumberOfRuns],dtype=tuple )
BigVirusArray = numpy.empty([NumberOfRuns],dtype=tuple )

ReallyBigVirusArray = numpy.empty([len(regenParams),NumberOfRuns],dtype=tuple )
BigMedianVirusArray = numpy.empty([len(regenParams)],dtype=tuple )

LongestTime = 0.0

for j in range(len(regenParams)):
    OUTERPATH = os.path.join(OuterPath).replace('\\','//')
    directory = os.path.join("".join(OUTERPATH+ "_" +str(regenParams[j])+ "_" +"Analysis"))
    if not os.path.exists(directory):
        os.makedirs(directory)
    OUTER_Path_to_Folder = os.path.abspath(directory)

    directory = os.path.join("".join(OUTERPATH+ "_" +"Graphs"))
    if not os.path.exists(directory):
        os.makedirs(directory)
    OUTEROUTER_Path_to_Folder = os.path.abspath(directory)

    for i in range(NumberOfRuns):
        os.system('cls' if os.name=='nt' else 'clear')
        print(str(regenParams[j])+"\n")
        print(i)
        InnerPath = r"/607_"+str(i)+"-FREECELL_-2.0-MOI_"+str(regenParams[j])+"-RP"


#        InnerPath = r"/607_"+str(i)+"-FREECELL_"+str(logMOI[j])+"-MOI_0.2-probi"
#        InnerPath = r"/7_"+str(i)+"-CELL2CELL_"+str(logMOI[j])+"-MOI"

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

        ReallyBigVirusArray[j][i] = VirusArray

        plt.rcParams['figure.figsize'] = [10, 10]
        fig, ax = plt.subplots()
        plt.plot(TimeArray, HealthyArray)
        plt.xticks(numpy.arange(0, TimeArray[-1], (TimeArray[-1]/10)))
        plt.tight_layout()
        plt.close()
        # fig.savefig(os.path.join(Inner_Path_to_Folder,"HealthyVsTime"+".pdf"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)
        fig.savefig(os.path.join(Inner_Path_to_Folder,"HealthyVsTime"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

        plt.rcParams['figure.figsize'] = [10, 10]
        fig, ax = plt.subplots()
        plt.plot(TimeArray, EclipseArray)
        plt.xticks(numpy.arange(0, TimeArray[-1], (TimeArray[-1]/10)))
        plt.tight_layout()
        plt.close()
        # fig.savefig(os.path.join(Inner_Path_to_Folder,"EclipseVsTime"+".p/df"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)
        fig.savefig(os.path.join(Inner_Path_to_Folder,"EclipseVsTime"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

        plt.rcParams['figure.figsize'] = [10, 10]
        fig, ax = plt.subplots()
        plt.plot(TimeArray, InfectedArray)
        plt.xticks(numpy.arange(0, TimeArray[-1], (TimeArray[-1]/10)))
        plt.tight_layout()
        plt.close()
        # fig.savefig(os.path.join(Inner_Path_to_Folder,"InfectedVsTime"+".pdf"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)
        fig.savefig(os.path.join(Inner_Path_to_Folder,"InfectedVsTime"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

        plt.rcParams['figure.figsize'] = [10, 10]
        fig, ax = plt.subplots()
        plt.plot(TimeArray, DeadArray)
        plt.xticks(numpy.arange(0, TimeArray[-1], (TimeArray[-1]/10)))
        plt.close()
        # fig.savefig(os.path.join(Inner_Path_to_Folder,"DeadVsTime"+".pdf"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)
        fig.savefig(os.path.join(Inner_Path_to_Folder,"DeadVsTime"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

        plt.rcParams['figure.figsize'] = [10, 10]
        fig, ax = plt.subplots()
        plt.yscale("log")
        plt.plot(TimeArray, VirusArray)
        plt.xticks(numpy.arange(0, TimeArray[-1], (TimeArray[-1]/10)))
        plt.tight_layout()
        plt.title("Virus vs Time")
        plt.ylabel("Amount of Virus")
        plt.xlabel("Time-(hr)")
        plt.close()
        # fig.savefig(os.path.join(Inner_Path_to_Folder,"VirusVsTime"+".pdf"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)
        fig.savefig(os.path.join(Inner_Path_to_Folder,"VirusVsTime"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

        BeginingTime = 2
        IndexOfPeakTime = numpy.argmax(BigVirusArray[i])

        BigPeakVirus[j][i] = max(BigVirusArray[i])
        BigTimeOfPeak[j][i] = BigTimeArray[i][IndexOfPeakTime]
        BigUpSlope[j][i] = growthfit(BigTimeArray[i], BigVirusArray[i], IndexOfPeakTime)[0]
        BigDownSlope[j][i] = 0.0 #decayfit(BigTimeArray[i], BigEclipseArray[i], BigInfectedArray[i], BigDeadArray[i], BigVirusArray[i])[0]
        BigEndTime[j][i] = BigTimeArray[i][-1]
        BigAUC[j][i] = numpy.trapz(numpy.log(BigVirusArray[i][2:]))
        BigChronicLoad[j][i] = numpy.average(VirusArray[500:600])

    MedianTimeArray = []
    MedianHealthyArray = []
    MedianEclipseArray = []
    MedianInfectedArray = []
    MedianDeadArray = []
    MedianVirusArray = []

    InbetweenTimeArray = []
    InbetweenHealthyArray = []
    InbetweenEclipseArray = []
    InbetweenInfectedArray = []
    InbetweenDeadArray = []
    InbetweenVirusArray = []

    for i in range(LongestTime):
        for k in range(NumberOfRuns):
            try:
                InbetweenTimeArray.append(BigTimeArray[k][i])
            except IndexError:
                InbetweenTimeArray.append(0.0)

            try:
                InbetweenHealthyArray.append(BigHealthyArray[k][i])
            except IndexError:
                InbetweenHealthyArray.append(0.0)

            try:
                InbetweenEclipseArray.append(BigEclipseArray[k][i])
            except IndexError:
                InbetweenEclipseArray.append(0.0)

            try:
                InbetweenInfectedArray.append(BigInfectedArray[k][i])
            except IndexError:
                InbetweenInfectedArray.append(0.0)

            try:
                InbetweenDeadArray.append(BigDeadArray[k][i])
            except IndexError:
                InbetweenDeadArray.append(0.0)

            try:
                InbetweenVirusArray.append(BigVirusArray[k][i])
            except IndexError:
                InbetweenVirusArray.append(0.0)

        MedianTimeArray.append(numpy.median(InbetweenTimeArray))
        MedianHealthyArray.append(numpy.median(InbetweenHealthyArray))
        MedianEclipseArray.append(numpy.median(InbetweenEclipseArray))
        MedianInfectedArray.append(numpy.median(InbetweenInfectedArray))
        MedianDeadArray.append(numpy.median(InbetweenDeadArray))
        MedianVirusArray.append(numpy.median(InbetweenVirusArray))

        InbetweenTimeArray.clear()
        InbetweenHealthyArray.clear()
        InbetweenEclipseArray.clear()
        InbetweenInfectedArray.clear()
        InbetweenDeadArray.clear()
        InbetweenVirusArray.clear()

    BigMedianVirusArray[j] = MedianVirusArray

    with open(os.path.join(OUTER_Path_to_Folder,str(regenParams[j])+"MedianPerTimeStep.txt"), 'w') as outfile:
        for i in range(len(MedianTimeArray)):
            print(str(int(MedianTimeArray[i]))+", "+str(int(MedianHealthyArray[i]))+", "+str(int(MedianEclipseArray[i]))+", "+str(int(MedianInfectedArray[i]))+", "+str(int(MedianDeadArray[i]))+", "+str(MedianVirusArray[i])+",",file=outfile)

    plt.rcParams['figure.figsize'] = [10, 10]
    fig, ax = plt.subplots()
    for i in BigHealthyArray:
        plt.plot(i, color = "gray", alpha = 0.3)
    plt.plot(MedianHealthyArray, color = "red")
    plt.xticks(numpy.arange(0, LongestTime, LongestTime/10))
    plt.tight_layout()
    plt.close()
    # fig.savefig(os.path.join(OUTER_Path_to_Folder,"regenParams-" + str(regenParams[j]) + "_Median_HealthyVsTime"+".pdf"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)
    fig.savefig(os.path.join(OUTER_Path_to_Folder,"regenParams-" + str(regenParams[j]) + "Median_HealthyVsTime"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

    plt.rcParams['figure.figsize'] = [10, 10]
    fig, ax = plt.subplots()
    for i in BigEclipseArray:
        plt.plot(i, color = "gray", alpha = 0.3)
    plt.plot(MedianEclipseArray, color = "red")
    plt.xticks(numpy.arange(0, LongestTime, LongestTime/10))
    plt.tight_layout()
    plt.close()
    # fig.savefig(os.path.join(OUTER_Path_to_Folder,"regenParams-" + str(regenParams[j]) + "Median_EclipseVsTime"+".pdf"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)
    fig.savefig(os.path.join(OUTER_Path_to_Folder,"regenParams-" + str(regenParams[j]) + "Median_EclipseVsTime"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

    plt.rcParams['figure.figsize'] = [10, 10]
    fig, ax = plt.subplots()
    for i in BigInfectedArray:
        plt.plot(i, color = "gray", alpha = 0.3)
    plt.plot(MedianInfectedArray, color = "red")
    plt.xticks(numpy.arange(0, LongestTime, LongestTime/10))
    plt.tight_layout()
    plt.close()
    # fig.savefig(os.path.join(OUTER_Path_to_Folder,"regenParams-" + str(regenParams[j]) + "Median_InfectedVsTime"+".pdf"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)
    fig.savefig(os.path.join(OUTER_Path_to_Folder,"regenParams-" + str(regenParams[j]) + "Median_InfectedVsTime"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

    plt.rcParams['figure.figsize'] = [10, 10]
    fig, ax = plt.subplots()
    for i in BigDeadArray:
        plt.plot(i, color = "gray", alpha = 0.3)
    plt.plot(MedianDeadArray, color = "red")
    plt.xticks(numpy.arange(0, LongestTime, LongestTime/10))
    plt.tight_layout()
    plt.close()
    # fig.savefig(os.path.join(OUTER_Path_to_Folder,"regenParams-" + str(regenParams[j]) + "Median_DeadVsTime"+".pdf"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)
    fig.savefig(os.path.join(OUTER_Path_to_Folder,"regenParams-" + str(regenParams[j]) + "Median_DeadVsTime"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

    plt.rcParams['figure.figsize'] = [10, 10]
    fig, ax = plt.subplots()
    plt.yscale("log")
    for i in BigVirusArray:
        plt.plot(i, color = "gray", alpha = 0.3)
    plt.plot(MedianVirusArray, color = "red")
    plt.xticks(numpy.arange(0, LongestTime, LongestTime/10))
    plt.tight_layout()
    plt.title("Median Virus vs Time")
    plt.ylabel("Median_ Virus Titer")
    plt.xlabel("Time-(hr)")
    plt.close()
    # fig.savefig(os.path.join(OUTER_Path_to_Folder,"regenParams-" + str(regenParams[j]) + "Median_VirusVsTime"+".pdf"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)
    fig.savefig(os.path.join(OUTER_Path_to_Folder,"regenParams-" + str(regenParams[j]) + "Median_VirusVsTime"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

#BigPeakVirus
    histbins = numpy.linspace(min(BigPeakVirus[j]),max(BigPeakVirus[j]),20)
    hist, bins = numpy.histogram(BigPeakVirus[j], bins=histbins, density=True)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2

    mean, STD = norm.fit(BigPeakVirus[j])
    x = numpy.linspace(min(BigPeakVirus[j]), max(BigPeakVirus[j]), 1000)
    p = norm.pdf(x, mean, STD)

    plt.rcParams['figure.figsize'] = [10, 10]
    fig, ax = plt.subplots()
    ax.bar(center, hist, align='center', width=width)
    plt.plot(x, p)
    plt.tight_layout()
    plt.close()
    # fig.savefig(os.path.join(OUTER_Path_to_Folder,str(regenParams[j])+"_PeakVirusHist"+".pdf"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)
    fig.savefig(os.path.join(OUTER_Path_to_Folder,str(regenParams[j])+"_PeakVirusHist"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

#BigTimeOfPeak
    histbins = numpy.linspace(min(BigTimeOfPeak[j]),max(BigTimeOfPeak[j]),20)
    hist, bins = numpy.histogram(BigTimeOfPeak[j], bins=histbins, density=True)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2

    mean, STD = norm.fit(BigTimeOfPeak[j])
    x = numpy.linspace(min(BigTimeOfPeak[j]), max(BigTimeOfPeak[j]), 1000)
    p = norm.pdf(x, mean, STD)

    plt.rcParams['figure.figsize'] = [10, 10]
    fig, ax = plt.subplots()
    ax.bar(center, hist, align='center', width=width)
    plt.plot(x, p)
    plt.ylabel("Frequency")
    plt.xlabel("Time of peak virus (hours)")
    plt.tight_layout()
    plt.close()
    # fig.savefig(os.path.join(OUTER_Path_to_Folder,str(regenParams[j])+"_TimeOfPeakHist"+".pdf"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)
    fig.savefig(os.path.join(OUTER_Path_to_Folder,str(regenParams[j])+"_TimeOfPeakHist"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

#BigUpSlope
    histbins = numpy.linspace(min(BigUpSlope[j]),max(BigUpSlope[j]),20)
    hist, bins = numpy.histogram(BigUpSlope[j], bins=histbins, density=True)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2

    mean, STD = norm.fit(BigUpSlope[j])
    x = numpy.linspace(min(BigUpSlope[j]), max(BigUpSlope[j]), 1000)
    p = norm.pdf(x, mean, STD)

    plt.rcParams['figure.figsize'] = [10, 10]
    fig, ax = plt.subplots()
    ax.bar(center, hist, align='center', width=width)
    plt.plot(x, p)
    plt.tight_layout()
    plt.close()
    # fig.savefig(os.path.join(OUTER_Path_to_Folder,str(regenParams[j])+"_UpSlopeHist"+".pdf"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)
    fig.savefig(os.path.join(OUTER_Path_to_Folder,str(regenParams[j])+"_UpSlopeHist"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

#BigDownSlope
    histbins = numpy.linspace(min(BigDownSlope[j]),max(BigDownSlope[j]),20)
    hist, bins = numpy.histogram(BigDownSlope[j], bins=histbins, density=True)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2

    mean, STD = norm.fit(BigDownSlope[j])
    x = numpy.linspace(min(BigDownSlope[j]), max(BigDownSlope[j]), 1000)
    p = norm.pdf(x, mean, STD)

    plt.rcParams['figure.figsize'] = [10, 10]
    fig, ax = plt.subplots()
    ax.bar(center, hist, align='center', width=width)
    plt.plot(x, p)
    plt.tight_layout()
    plt.close()
    # fig.savefig(os.path.join(OUTER_Path_to_Folder,str(regenParams[j])+"_DownSlopeHist"+".pdf"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)
    fig.savefig(os.path.join(OUTER_Path_to_Folder,str(regenParams[j])+"_DownSlopeHist"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

#BigEndTime
    histbins = numpy.linspace(min(BigEndTime[j]),max(BigEndTime[j]),20)
    hist, bins = numpy.histogram(BigEndTime[j], bins=histbins, density=True)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2

    mean, STD = norm.fit(BigEndTime[j])
    x = numpy.linspace(min(BigEndTime[j]), max(BigEndTime[j]), 1000)
    p = norm.pdf(x, mean, STD)

    plt.rcParams['figure.figsize'] = [10, 10]
    fig, ax = plt.subplots()
    ax.bar(center, hist, align='center', width=width)
    plt.plot(x, p)
    plt.tight_layout()
    plt.close()
    # fig.savefig(os.path.join(OUTER_Path_to_Folder,str(regenParams[j])+"_DownSlopeHist"+".pdf"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)
    fig.savefig(os.path.join(OUTER_Path_to_Folder,str(regenParams[j])+"_DownSlopeHist"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

#BigAUC
    histbins = numpy.linspace(min(BigAUC[j]),max(BigAUC[j]),20)
    hist, bins = numpy.histogram(BigAUC[j], bins=histbins, density=True)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2

    mean, STD = norm.fit(BigAUC[j])
    x = numpy.linspace(min(BigAUC[j]), max(BigAUC[j]), 1000)
    p = norm.pdf(x, mean, STD)

    plt.rcParams['figure.figsize'] = [10, 10]
    fig, ax = plt.subplots()
    ax.bar(center, hist, align='center', width=width)
    plt.plot(x, p)
    plt.tight_layout()
    plt.close()
    # fig.savefig(os.path.join(OUTER_Path_to_Folder,str(regenParams[j])+"_DownSlopeHist"+".pdf"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)
    fig.savefig(os.path.join(OUTER_Path_to_Folder,str(regenParams[j])+"_DownSlopeHist"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

#BigChronicLoad
    histbins = numpy.linspace(min(BigChronicLoad[j]),max(BigChronicLoad[j]),20)
    hist, bins = numpy.histogram(BigChronicLoad[j], bins=histbins, density=True)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2

    mean, STD = norm.fit(BigChronicLoad[j])
    x = numpy.linspace(min(BigChronicLoad[j]), max(BigChronicLoad[j]), 1000)
    p = norm.pdf(x, mean, STD)

    plt.rcParams['figure.figsize'] = [10, 10]
    fig, ax = plt.subplots()
    ax.bar(center, hist, align='center', width=width)
    plt.plot(x, p)
    plt.tight_layout()
    plt.close()
    # fig.savefig(os.path.join(OUTER_Path_to_Folder,str(regenParams[j])+"_DownSlopeHist"+".pdf"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)
    fig.savefig(os.path.join(OUTER_Path_to_Folder,str(regenParams[j])+"_DownSlopeHist"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

#Make Aspect vs regenParams Graphs
    size = 32
    capsizes = [50, 21, 50, 21, 50, 21, 50, 21, 50, 21]
    axessize = 2

    plt.rcParams['xtick.labelsize'] = size
    plt.rcParams['ytick.labelsize'] = size
    plt.rcParams['axes.labelsize'] = size
    plt.rcParams['figure.figsize'] = [12 , 8]
    plt.rcParams['xtick.major.width'] = axessize
    plt.rcParams['xtick.minor.width'] = axessize
    plt.rcParams['ytick.major.width'] = axessize
    plt.rcParams['ytick.minor.width'] = axessize
    plt.rc('axes', linewidth=axessize)

    fig, ax = plt.subplots()
    ax.yaxis.set_major_locator(ticker.MaxNLocator(6))
    holder = []
    holder1 = []
    for i in range(len(regenParams)):
        holder.append(numpy.mean(BigPeakVirus[i][:]))
        holder1.append(numpy.std(BigPeakVirus[i][:]))
    plt.errorbar(regenParamsFloat, holder, yerr=holder1, fmt = "-o", ecolor="black", capsize=10, linewidth=axessize)
    plt.xlabel("Regeneration Rate")
    plt.ylabel("Peak Virus")
    plt.xscale("log")
    plt.yscale("log")
    plt.tight_layout()
    plt.close()
    # fig.savefig(os.path.join(OUTEROUTER_Path_to_Folder, "PeakViralTitter"+".pdf"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)
    fig.savefig(os.path.join(OUTEROUTER_Path_to_Folder, "PeakViralTitter"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

    fig, ax = plt.subplots()
    ax.yaxis.set_major_locator(ticker.MaxNLocator(6))
    holder = []
    holder1 = []
    for i in range(len(regenParams)):
        holder.append(numpy.mean(BigTimeOfPeak[i][:])/24)
        holder1.append(numpy.std(BigTimeOfPeak[i][:])/24)
    plt.errorbar(regenParamsFloat, holder, yerr=holder1, fmt = "-o", ecolor="black", capsize=10)
    plt.xlabel("Regeneration Rate")
    plt.ylabel("Time of Peak (day)")
    plt.xscale("log")
    plt.tight_layout()
    plt.close()
    # fig.savefig(os.path.join(OUTEROUTER_Path_to_Folder, "TimeofPeakViralTitter"+".pdf"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)
    fig.savefig(os.path.join(OUTEROUTER_Path_to_Folder, "TimeofPeakViralTitter"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

    fig, ax = plt.subplots()
    ax.yaxis.set_major_locator(ticker.MaxNLocator(6))
    holder = []
    holder1 = []
    for i in range(len(regenParams)):
        holder.append(numpy.mean(BigUpSlope[i][:])*24)
        holder1.append(numpy.std(BigUpSlope[i][:])*24)
    plt.errorbar(regenParamsFloat, holder, yerr=holder1, fmt = "-o", ecolor="black", capsize=10)
    plt.xlabel("Regeneration Rate")
    plt.ylabel("Growth Rate (virus/day)")
    plt.xscale("log")
    plt.tight_layout()
    plt.close()
    # fig.savefig(os.path.join(OUTEROUTER_Path_to_Folder, "UpSlope"+".pdf"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)
    fig.savefig(os.path.join(OUTEROUTER_Path_to_Folder, "UpSlope"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

    # fig, ax = plt.subplots()
    # ax.yaxis.set_major_locator(ticker.MaxNLocator(6))
    # holder = []
    # holder1 = []
    # holder2 = []
    # holder3 = []
    # for i in range(len(regenParams)):
    #     holder.append(numpy.mean(BigDownSlope[i][:])*24)
    #     holder1.append(numpy.std(BigDownSlope[i][:])*24)
    #     holder2.append(numpy.round(numpy.mean(BigDownSlope[i][:]),3)*24)
    #     holder3.append(numpy.round(numpy.std(BigDownSlope[i][:]),3)*24)
    # plt.errorbar(regenParamsFloat, holder, yerr=holder1, fmt = "-o", ecolor="black", capsize=10)
    # plt.xlabel("Regeneration Rate")
    # plt.ylabel("Decay Rate (virus/day)")
    # plt.tight_layout()
    # plt.close()
    # fig.savefig(os.path.join(OUTEROUTER_Path_to_Folder, "DownSlope"+".pdf"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)
    # fig, ax = plt.subplots()
    # ax.yaxis.set_major_locator(ticker.MaxNLocator(6))
    # plt.errorbar(regenParamsFloat, holder2, yerr=holder3, fmt = "-o", ecolor="black", capsize=10)
    # plt.xlabel("Regeneration Rate")
    # plt.ylabel("Decay Rate (virus/day)")
    # plt.tight_layout()
    # plt.close()
    # fig.savefig(os.path.join(OUTEROUTER_Path_to_Folder, "DownSlopeRound"+".pdf"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)
    # fig.savefig(os.path.join(OUTEROUTER_Path_to_Folder, "DownSlopeRound"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

    # fig, ax = plt.subplots()
    # ax.yaxis.set_major_locator(ticker.MaxNLocator(6))
    # holder = []
    # holder1 = []
    # for i in range(len(regenParams)):
    #     holder.append(numpy.mean(BigEndTime[i][:])/24)
    #     holder1.append(numpy.std(BigEndTime[i][:])/24)
    # plt.errorbar(regenParamsFloat, holder, yerr=holder1, fmt = "-o", ecolor="black", capsize=10)
    # plt.xlabel("Regeneration Rate")
    # plt.ylabel("Infection Duration (day)")
    # plt.xscale("log")
    # plt.tight_layout()
    # plt.close()
    # fig.savefig(os.path.join(OUTEROUTER_Path_to_Folder, "InfectionDuration"+".pdf"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)
    # fig.savefig(os.path.join(OUTEROUTER_Path_to_Folder, "InfectionDuration"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

    fig, ax = plt.subplots()
    ax.yaxis.set_major_locator(ticker.MaxNLocator(6))
    holder = []
    holder1 = []
    for i in range(len(regenParams)):
        holder.append(numpy.mean(BigAUC[i][:])/24)
        holder1.append(numpy.std(BigAUC[i][:])/24)
    plt.errorbar(regenParamsFloat, holder, yerr=holder1, fmt = "-o", ecolor="black", capsize=10)
    plt.xlabel("Regeneration Rate")
    plt.ylabel("AUC (log(virus)*day)")
    plt.xscale("log")
    plt.tight_layout()
    plt.close()
    # fig.savefig(os.path.join(OUTEROUTER_Path_to_Folder, "AUC"+".pdf"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)
    fig.savefig(os.path.join(OUTEROUTER_Path_to_Folder, "AUC"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

    fig, ax = plt.subplots()
    ax.yaxis.set_major_locator(ticker.MaxNLocator(6))
    holder = []
    holder1 = []
    for i in range(len(regenParams)):
        holder.append(numpy.mean(BigChronicLoad[i][:]))
        holder1.append(numpy.std(BigChronicLoad[i][:]))
    plt.errorbar(regenParamsFloat, holder, yerr=holder1, fmt = "-o", ecolor="black", capsize=10)
    plt.xlabel("Regeneration Rate")
    plt.ylabel("Chronic Load")
    plt.xscale("log")
    plt.tight_layout()
    plt.close()
    # fig.savefig(os.path.join(OUTEROUTER_Path_to_Folder, "ChronicLoad"+".pdf"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)
    fig.savefig(os.path.join(OUTEROUTER_Path_to_Folder, "ChronicLoad"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

    with open(os.path.join(OUTEROUTER_Path_to_Folder,"PeakVirus.txt"), 'w') as outfile:
        for i in range(len(regenParams)):
            print(str(numpy.mean(BigPeakVirus[i][:]))+","+str(numpy.std(BigPeakVirus[i][:])),file=outfile)

    with open(os.path.join(OUTEROUTER_Path_to_Folder,"TimeOfPeakVirus.txt"), 'w') as outfile:
        for i in range(len(regenParams)):
            print(str(numpy.mean(BigTimeOfPeak[i][:]))+","+str(numpy.std(BigTimeOfPeak[i][:])),file=outfile)

    with open(os.path.join(OUTEROUTER_Path_to_Folder,"UpSlope.txt"), 'w') as outfile:
        for i in range(len(regenParams)):
            print(str(numpy.mean(BigUpSlope[i][:]))+","+str(numpy.std(BigUpSlope[i][:])),file=outfile)
    # with open(os.path.join(OUTEROUTER_Path_to_Folder,"DownSlope.txt"), 'w') as outfile:
    #     for i in range(len(regenParams)):
    #         print(str(numpy.mean(BigDownSlope[i][:]))+","+str(numpy.std(BigDownSlope[i][:])),file=outfile)
    # with open(os.path.join(OUTEROUTER_Path_to_Folder,"InfectionDuration.txt"), 'w') as outfile:
    #     for i in range(len(regenParams)):
    #         print(str(numpy.mean(BigEndTime[i][:]))+","+str(numpy.std(BigEndTime[i][:])),file=outfile)
    with open(os.path.join(OUTEROUTER_Path_to_Folder,"AUC.txt"), 'w') as outfile:
        for i in range(len(regenParams)):
            print(str(numpy.mean(BigAUC[i][:]))+","+str(numpy.std(BigAUC[i][:])),file=outfile)
    with open(os.path.join(OUTEROUTER_Path_to_Folder,"AUC.txt"), 'w') as outfile:
        for i in range(len(regenParams)):
            print(str(numpy.mean(BigChronicLoad[i][:]))+","+str(numpy.std(BigChronicLoad[i][:])),file=outfile)

    BigArray = [BigPeakVirus, BigTimeOfPeak, BigUpSlope, BigDownSlope, BigEndTime, BigAUC, BigChronicLoad]
    BigStrArray = ["BigPeakVirus", "BigTimeOfPeak", "BigUpSlope", "BigDownSlope", "BigEndTime", "BigAUC", "BigChronicLoad"]

    for q in range(len(BigArray)):
        #    Shapiro-Wilk
        SStat, SP = shapiro(BigArray[q][j])
        #    D'Agostino's
        DStat, DP = normaltest(BigArray[q][j])
        #    Anderson-Darling
        result = anderson(BigArray[q][j])

        with open(os.path.join(OUTER_Path_to_Folder,str(regenParams[j])+"_Normality"+".txt"),'a') as outfile:

            print(str(BigStrArray[q]),file=outfile)

        #    Shapiro-Wilk
            print("Shapiro-Wilk",file=outfile)
            print("Statistics= "+"{0:.3g}".format(SStat)+", p= "+"{0:.3g}".format(SP),file=outfile)
            SAlpha = 0.05
            if SP > SAlpha:
	            print('Sample looks Gaussian (fail to reject H0)',file=outfile)
            else:
	            print('Sample does not look Gaussian (reject H0)',file=outfile)

        #    D'Agostino's
            print("D'Agostino's",file=outfile)
            print("Statistics= "+"{0:.3g}".format(DStat)+", p= "+"{0:.3g}".format(DP),file=outfile)
            DAlpha = 0.05
            if DP > DAlpha:
	            print('Sample looks Gaussian (fail to reject H0)',file=outfile)
            else:
	            print('Sample does not look Gaussian (reject H0)',file=outfile)

        #    Anderson-Darling
            print("Anderson-Darling",file=outfile)
            print("Statistic: "+"{0:.3g}".format(result.statistic),file=outfile)
            p = 0
            for i in range(len(result.critical_values)):
	            sl, cv = result.significance_level[i], result.critical_values[i]
	            if result.statistic < result.critical_values[i]:
		            print("{0:.3g}".format(sl)+":"+"{0:.3g}".format(cv)+", data looks normal (fail to reject H0)",file=outfile)
	            else:
		            print("{0:.3g}".format(sl)+":"+"{0:.3g}".format(cv)+", data does not looks normal (reject H0)",file=outfile)

            print("\n",file=outfile)





#directory = os.path.join("".join(OUTERPATH+ "_" +"AllOnOne"))
#if not os.path.exists(directory):
#    os.makedirs(directory)
#AllOnOne_Path_to_Folder = os.path.abspath(directory)
##All on one
#size = 32
#plt.rcParams['xtick.labelsize'] = size
#plt.rcParams['ytick.labelsize'] = size
#plt.rcParams['axes.labelsize'] = size
#plt.rc('legend', fontsize=size*(7/10))
#plt.rc('legend', title_fontsize=size*(7/10))
#plt.rcParams['figure.figsize'] = [15, 15]
#fig, ax = plt.subplots()
#plt.yscale("log")
##plt.ylim(numpy.power(10,5), numpy.power(10,9))
##plt.xlim(0, 60*24)
##colors = ["Purple", "Purple", "Orange", "Orange", "Blue", "Blue", "Green", "Green", "Red", "Red"]
#colors = ["Cyan", "Purple", "Orange", "Blue", "Green", "Red"]
##alphas = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.05, 0.05]
#alphas = [0.1, 0.1, 0.1, 0.1, 0.1, 0.05]
#for j in range(len(MOI)):
##    if(j%2==1):
#    for i in ReallyBigVirusArray[j]:
#        plt.plot(i, color = colors[j], alpha = alphas[j])
#    plt.plot(BigMedianVirusArray[j], color = colors[j], label = MOILabels[j])
##plt.xticks(numpy.linspace(0, 60*24, 16), numpy.linspace(0, 60*24, 16, dtype = int)//24)
#plt.tight_layout()
##plt.title("Median Virus vs Time")
#plt.ylabel("Median Virus Titer")
#plt.xlabel("Time (day)")
#plt.legend(title="MOI")
#plt.close()
#fig.savefig(os.path.join(AllOnOne_Path_to_Folder,"AllOnOne_Median_VirusVsTime"+".pdf"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)
#fig.savefig(os.path.join(AllOnOne_Path_to_Folder,"AllOnOne_Median_VirusVsTime"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

#print("Data Analysis Complete")
