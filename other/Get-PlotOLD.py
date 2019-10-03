import numpy
import datetime
import matplotlib.pyplot as plt
import os

os.system('cls' if os.name=='nt' else 'clear')

s = r"C:\Users\BaylorFain\Documents\Python\Research\Viral_Model\ViralModel\FREECELL\-4.0-MOI\(4-FREECELL)(-4.0-MOI)(2019-04-05)(12-44-11)"
INITIALPATH = os.path.join(s).replace('\\','//')

PRINT = 0
ANALYSIS = 1

NumberOfLayers = 193 #607 is a million hexigon in a circle
Nt = int(376)

if ANALYSIS == 1:
    directory = os.path.join("".join(INITIALPATH + "/Analysis"))
    if not os.path.exists(directory):
        os.makedirs(directory)
    Path_to_Folder = os.path.abspath(directory) # Figures out the absolute path for you in case your working directory moves around. 
    
    Nx=(2*NumberOfLayers)-1
    Ny=(2*NumberOfLayers)-1
    N=Nx*Ny
    
    HealthyArray = numpy.empty([Nt,1],dtype=float )
    EclipseArray = numpy.empty([Nt,1],dtype=float )
    InfectedArray = numpy.empty([Nt,1],dtype=float )
    DeadArray = numpy.empty([Nt,1],dtype=float )
    VirusArray = numpy.empty([Nt,1],dtype=float )

    HealthyArray.fill(0)
    EclipseArray.fill(0)
    InfectedArray.fill(0)
    DeadArray.fill(0)
    VirusArray.fill(0)
    
    with open(os.path.join(INITIALPATH,"cells_over_time.txt")) as infile:
        CellData = infile.read()
        CellData = CellData.replace("\n",",")
        CellData = CellData.replace(" ","")
        CellData = CellData.split(",")
    with open(os.path.join(INITIALPATH,"virus_over_time.txt")) as infile:
        VirusData = infile.read()
        VirusData = VirusData.replace("\n","")
        VirusData = VirusData.replace(" ","")
        VirusData = VirusData.split(",")
    
    for q in range(Nt):
        HealthCount = 0
        Eclipsecount = 0
        InfectedCount = 0
        DeadCount = 0
        VirusCount = 0
        for n in range(Nx):
            for i in range(Nx):
                v = VirusData[(i*Nx)+n+N*q]
                VirusCount += float(v)
                c = CellData[(i*Nx)+n+N*q]
                #if c != "0":
                if c == "h":
                    HealthCount += 1
                elif c == "e":
                    Eclipsecount += 1
                elif c == "i":
                    InfectedCount += 1
                elif c == "d":
                    DeadCount += 1
                elif c == "0":
                    "Nothing"
        HealthyArray[q] = HealthCount
        EclipseArray[q] = Eclipsecount
        InfectedArray[q] = InfectedCount
        DeadArray[q] = DeadCount
        VirusArray[q] = VirusCount
    
    with open(os.path.join(Path_to_Folder,"healthy_cell_analysis.txt"),'a') as outfile:
        for RowIndex in HealthyArray:
            print("".join(["{}".format(ColumnIndex) for ColumnIndex in RowIndex]),
                  file=outfile)
    with open(os.path.join(Path_to_Folder,"eclipse_cell_analysis.txt"),'a') as outfile:
        for RowIndex in EclipseArray:
            print("".join(["{}".format(ColumnIndex) for ColumnIndex in RowIndex]),
                  file=outfile)
    with open(os.path.join(Path_to_Folder,"infected_cell_analysis.txt"),'a') as outfile:
        for RowIndex in InfectedArray:
            print("".join(["{}".format(ColumnIndex) for ColumnIndex in RowIndex]),
                  file=outfile)
    with open(os.path.join(Path_to_Folder,"dead_cell_analysis.txt"),'a') as outfile:
        for RowIndex in DeadArray:
            print("".join(["{}".format(ColumnIndex) for ColumnIndex in RowIndex]),
                  file=outfile)
    with open(os.path.join(Path_to_Folder,"virus_analysis.txt"),'a') as outfile:
        for RowIndex in VirusArray:
            print("".join(["{}".format(ColumnIndex) for ColumnIndex in RowIndex]),
                  file=outfile)
      
    plt.rcParams['figure.figsize'] = [10, 10]
    fig, ax = plt.subplots()
    plt.plot(HealthyArray)
    plt.close()
    fig.savefig(os.path.join(Path_to_Folder,"HealthyVsTime"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)
    
    plt.rcParams['figure.figsize'] = [10, 10]
    fig, ax = plt.subplots()
    plt.plot(EclipseArray)
    plt.close()
    fig.savefig(os.path.join(Path_to_Folder,"EclipseVsTime"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)
    
    plt.rcParams['figure.figsize'] = [10, 10]
    fig, ax = plt.subplots()
    plt.plot(InfectedArray)
    plt.close()
    fig.savefig(os.path.join(Path_to_Folder,"InfectedVsTime"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)
    
    plt.rcParams['figure.figsize'] = [10, 10]
    fig, ax = plt.subplots()
    plt.plot(DeadArray)
    plt.close()
    fig.savefig(os.path.join(Path_to_Folder,"DeadVsTime"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)
    
    plt.rcParams['figure.figsize'] = [10, 10]
    fig, ax = plt.subplots()
    plt.yscale("log")
    plt.plot(VirusArray)
    plt.title("Virus vs Time")
    plt.ylabel("Amount of Virus")
    plt.xlabel("Time-(hr)")
    plt.close()
    fig.savefig(os.path.join(Path_to_Folder,"VirusVsTime"+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

    print("Data Analysis Complete")
    os.startfile(Path_to_Folder)
    
if PRINT == 1:
    s = (2.0/3.0)
    RadiusScale = 0
    for i in range(NumberOfLayers):
        if i == 0:
            RadiusScale = RadiusScale + 1
        else:
            if (i)%2 == 1:
                RadiusScale = RadiusScale + 1
            else:
                RadiusScale = RadiusScale + 2
    RadiusOfCircle = s*RadiusScale

    count = 0
    for i in range(NumberOfLayers):
        count = count + i
    N=(count)*6+1

    coord = numpy.empty([N,3],dtype=float )
    percyclecoord = numpy.empty([N,3],dtype=tuple)

    percyclecoord.fill(0)
    coord.fill(0)

    percyclecoord[0] = [0,0,0]
    for j in range(NumberOfLayers):
        for i in range(2*j):
            if i < j:
                temp = i
            percyclecoord[i+(j-1)*j+1][0] =  -temp-1
            percyclecoord[i+(j-1)*j+1][1] =   temp+j-i
            percyclecoord[i+(j-1)*j+1][2] =  -j+1+i

    c0 = [percyclecoord[0][0], percyclecoord[0][1], percyclecoord[0][2]]
    coord[0][2] = c0[2]
    coord[0][1] = c0[1]
    coord[0][0] = c0[0]

    count = 0
    for j in range(int(N/3)):
        for i in range(3):
            coord[(i+0)%3+3*j+1][2] = percyclecoord[j+1][i]+c0[i]
            coord[(i+1)%3+3*j+1][1] = percyclecoord[j+1][i]+c0[i]
            coord[(i+2)%3+3*j+1][0] = percyclecoord[j+1][i]+c0[i]

    hi = coord[0][0]
    vi = coord[0][2]
    xmin = numpy.Inf
    for i in range(len(coord)):
        xcoord = coord[i][0]
        if coord[i][0] < xmin:
            xmin = coord[i][0]
        ycoord = (2.0 * numpy.sin(numpy.radians(60)) * (coord[i][1] - coord[i][2]) /3.0)+vi
        dist = numpy.sqrt((xcoord-hi)**2+(ycoord-vi)**2)
        if dist >= RadiusOfCircle:
            coord[i][0] = 5000.0
            coord[i][1] = 0.0
            coord[i][2] = 0.0

    QRIndexing = numpy.empty([(2*NumberOfLayers)-1,(2*NumberOfLayers)-1],dtype=float)
    QRIndexing.fill(0.0)
    for i in range(len(coord)):
        if coord[i][0] != 5000:
            QRIndexing[int(coord[i][2])-int(xmin)][int(coord[i][0])-int(xmin)] = 1.0
    Index = numpy.where(QRIndexing != 0)
    CoordanateNested = numpy.empty([(2*NumberOfLayers)-1,(2*NumberOfLayers)-1],dtype=tuple)
    CoordanateNested.fill(0)
    for i in range(len(Index[0])):
        CoordanateNested[Index[0][i]][Index[1][i]] = [Index[0][i]+xmin,-(Index[1][i]+xmin+Index[0][i]+xmin),Index[1][i]+xmin]
        
        
    print("Printing Biginnings")
        
    Nx=(2*NumberOfLayers)-1
    Ny=(2*NumberOfLayers)-1
    N=Nx*Ny

    directory = os.path.join("".join(INITIALPATH + "/Photos"))
    if not os.path.exists(directory):
        os.makedirs(directory)
    Path_to_Folder = os.path.abspath(directory) # Figures out the absolute path for you in case your working directory moves around. 

    hcoord = numpy.empty([Nt,((2*NumberOfLayers)-1)*((2*NumberOfLayers)-1)],dtype=float )
    vcoord = numpy.empty([Nt,((2*NumberOfLayers)-1)*((2*NumberOfLayers)-1)],dtype=float )
    colors = numpy.empty([Nt,((2*NumberOfLayers)-1)*((2*NumberOfLayers)-1)],dtype=tuple )
    virus = numpy.empty([Nt,((2*NumberOfLayers)-1)*((2*NumberOfLayers)-1)],dtype=float )
    hcoord.fill(float(NumberOfLayers))
    vcoord.fill(float(NumberOfLayers))
    colors.fill((1.0, 1.0, 1.0))
    virus.fill(0.0)

    with open(os.path.join(INITIALPATH,"cells_over_time.txt")) as infile:
        CellData = infile.read()
        CellData = CellData.replace("\n",",")
        CellData = CellData.replace(" ","")
        CellData = CellData.split(",")
    with open(os.path.join(INITIALPATH,"virus_over_time.txt")) as infile:
        VirusData = infile.read()
        VirusData = VirusData.replace("\n","")
        VirusData = VirusData.replace(" ","")
        VirusData = VirusData.split(",")

    maxdata1 = 0.0 
    mindata1 = 100000000.0
    for i in range(N*Nt):
        number = float(VirusData[i])
        if number>maxdata1:
            maxdata1 = number
        elif number<mindata1:
            mindata1 = number

    count = 0
    NumberOfPics = Nt
    for q in range(Nt):
        for n in range(Nx):
            for i in range(Nx):
                v = VirusData[(i*Nx)+n+N*q]
                virus[q][(i*Nx)+n] = (float(v)-mindata1)/(maxdata1-mindata1)
                
                c = CellData[(i*Nx)+n+N*q]
                #if c != "0":
                if c == "h":
                    color = (0.0, 1.0, 0.0)
                elif c == "e":
                    color = (0.0, 1.0, 1.0)
                elif c == "i":
                    color = (1.0, 0.0, 0.0)
                elif c == "d":
                    color = (0.0, 0.0, 0.0)
                elif c == "0":
                    color = (1.0, 1.0, 1.0)
                colors[q][(i*Nx)+n] = color
                
                if CoordanateNested[i][n] != 0:
                    # Horizontal cartesian coords
                    hcoord[q][(i*Nx)+n] = CoordanateNested[i][n][0]
                    # Vertical cartersian coords
                    vcoord[q][(i*Nx)+n] = (2.0 * numpy.sin(numpy.radians(60))*
                                           (CoordanateNested[i][n][1] - CoordanateNested[i][n][2]) /3.0)
                    
        plt.rcParams['figure.figsize'] = [10, 10]
        fig, ax = plt.subplots()
        plt.axis("off")
#         ax.scatter(hcoord[q], vcoord[q], facecolor=colors[q], edgecolors=colors[q], alpha=0.3, s = 1)
        ax.scatter(hcoord[q], vcoord[q], facecolor=colors[q], edgecolors=colors[q], alpha=0.3, s = 2250, marker = "H")
#         ax.scatter(hcoord[q], vcoord[q], facecolor=(0.5,0.0,0.5), edgecolors=(0.5,0.0,0.5),alpha=0.5, s = 2250*virus[q], marker = "H")
        if q == 0:
            [ymin, ymax] = ax.get_ylim()
            [xmin, xmax] = ax.get_xlim()
        ax.set_ylim(ymin, ymax)
        ax.set_xlim(xmin, xmax)
        plt.close()
        fig.savefig(os.path.join(Path_to_Folder,"cell"+"{}".format(q)+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)
        print(q)
        count = count+1
        if count >= NumberOfPics:
            break
    print("Image Printing Complete")
    os.startfile(Path_to_Folder)

if (ANALYSIS == 0) and (PRINT == 0):
    print("Nothing")
