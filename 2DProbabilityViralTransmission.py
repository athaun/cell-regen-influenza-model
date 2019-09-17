################################################################################
#
#                          Virus Model
#
################################################################################
#Ubuntu - clt + m to comment out mutiple lines
#Jupyter - clt + / to comment out mutiple lines
################################################################################
# Importing modules and packages
################################################################################
import numpy
import datetime
import gc
import sys
import os

################################################################################
#        Clear Terminal
################################################################################
os.system('cls' if os.name=='nt' else 'clear')

################################################################################
# Defining functions
################################################################################
def Te():
    """
    Picks a random number from the gamma distribution
    The number is to be used as a time step in the Eclipse Time Matrix
    """
    return numpy.random.gamma(TauE, TauE/numpy.sqrt(ne))

def Ti():
    """
    Picks a random number from the gamma distribution
    The number is to be used as a time step in the Infected Time Matrix
    """
    return numpy.random.gamma(TauI, TauI/numpy.sqrt(ni))

def PU1():
    """
    Picks a random number from a uniform distribution
    This probability
    """
    return numpy.random.uniform(low=0.0, high=1.0, size=None)

def exponentialDistro():
    """
    TODO: create exponential distrobution with a mean of 720 somehow uses to determine cell regen
    """
    return numpy.random.exponential(3 * 24)

# for printing different colored text to the terminal
class Tc:
    INFO = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    END = '\033[0m' # place at end of colored text to stop "color bleed"
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def printInfo(text):
    print(TextColors.INFO + text + TextColors.END)

def printWarning(text):
    print(TextColors.RED + TextColors.BOLD + TextColors.UNDERLINE + text + TextColors.END)


################################################################################
#      Simulation Parameters
################################################################################
CELL2CELL = 0
FREECELL = 1
timestep = 0.005    # Time step for model (No larger than 0.01 hour) 0.005 hr = 18 sec, (1/3600) hr = 1 sec
endtime = 365*24   # Days in hours, Need 7 days
Save = (1/timestep) #the number of time the program saves to file, (1/timestep) results in 1 save very simulated hour
NumberOfLayers = 7 #607 is a million hexigon in a circle
NumberOfRuns = 2

################################################################################
#    Physical Parameters
################################################################################
MOI = 10**(-2) #(10**-5) to 1 # fraction of cells intially infected
beta = 2.0 # infection rate
rho = 562800 # the probability of a cell getting inected (higher = more likely)
D = 6*10**(-12) # Diffusion coefficient (spaceUnit^2/seconds) how quickly the virus spread out
c = 0.105 # clearance rate of the virus (virus cleared/time unit)
deltx = 25.0e-06
deltxprime = deltx*2
Dtsx2 = D*timestep*(deltxprime**-2)

################################################################################
#    Probability Variables
################################################################################
TauI = 12.0   # Avg time for infection
TauE = 6.0    # Avg time for eclipse
ne = 30.0     #
ni = 100.0     #
probi = 0.2  # Probability per unit time of cell to cell infection (/hour)






totalRegenerations = 0






################################################################################

#                Loop For The number Of Simulations To Run Per Setting

################################################################################
for BigIndex in range(NumberOfRuns):
    ############################################################################
    #        Creating Save Path
    ############################################################################
    if CELL2CELL == 1:
        s = "CELL2CELL"
        if FREECELL == 1:
            s = "Both"
    elif  CELL2CELL == 0:
        s = "FREECELL"
        if FREECELL == 0:
            s = "Neither"
    directory = os.path.join("".join("ViralModel/"+"{}".format(NumberOfLayers)+"_{}".format(BigIndex)+"-"+s+"_{}-".format(numpy.log10(MOI))+"MOI"))

    if not os.path.exists(directory):
        os.makedirs(directory)
    Path_to_Folder = os.path.abspath(directory) # Figures out the absolute path for you in case your working directory moves around.

    ############################################################################
    #        Creating placeholder variables for multipy runs
    ############################################################################
    cell2cell = CELL2CELL
    freecell = FREECELL

    ############################################################################
    #          Building Cells
    ############################################################################
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

    QRIndexing = numpy.empty([(2*NumberOfLayers)-1,(2*NumberOfLayers)-1],dtype=str)
    QRIndexing.fill("0")
    for i in range(len(coord)):
        if coord[i][0] != 5000:
            QRIndexing[int(coord[i][2])-int(xmin)][int(coord[i][0])-int(xmin)] = "h"
    Index = numpy.where(QRIndexing != "0")
    CoordanateNested = numpy.empty([(2*NumberOfLayers)-1,(2*NumberOfLayers)-1],dtype=tuple)
    CoordanateNested.fill(0)
    for i in range(len(Index[0])):
        CoordanateNested[Index[0][i]][Index[1][i]] = [Index[0][i]+xmin,-(Index[1][i]+xmin+Index[0][i]+xmin),Index[1][i]+xmin]

    with open(os.path.join(Path_to_Folder,"InitialCellLocations.txt"),"w") as outfile:
        for n in QRIndexing:
            print(", ".join(["{}".format(i) for i in n]),file = outfile)

    with open(os.path.join(Path_to_Folder,"Parameters.txt"),"w") as outfile:
        # infoString = str(Tc.BLUE + "Day\t" + Tc.GREEN + "\t Healthy\t" + Tc.YELLOW + "Ecl/Inf\t\t" + Tc.RED + "Dead\t\t" + Tc.INFO + Tc.BOLD + "# Regens" + Tc.END)
        # print(infoString);

        print("Day \t",
              "NumberHealthy\t",
              "NumberEclipse\t",
              "NumberInfected\t",
              "NumberDead\t",
              "totalRegenerations")

        print("Hexagon Side Length =",s,
              "\nNumber of Layers =",NumberOfLayers,
              "\nRadius of Circle =",RadiusOfCircle,
              "\nNumber of Cells =", len(Index[0]),
              file = outfile)

    LocationData = QRIndexing
    NumberOfLayers = int(NumberOfLayers)
    NumberOfCells = float(len(Index[0]))

    ############################################################################
    #    Number of Cells
    ############################################################################
    Ni = NumberOfCells*MOI         #Number of infected cells
    if(Ni < 1): print("Use larger MOI"); sys.exit(0) #Exits program if initial number of infected cells is less then one.
    Nx = (2*NumberOfLayers-1)      # Range of cells on x axis
    Ny = (2*NumberOfLayers-1)      # Range of cells on y axis

    ############################################################################
    # Makeing empty matrices
    ############################################################################
    cells = numpy.empty([Ny,Nx],dtype = str)  #Produces a matrix for the cells
    ecl = numpy.empty([Ny,Nx],dtype = float) #Produces a time matrix for after eclipse phase (e)
    inf = numpy.empty([Ny,Nx],dtype = float) #Produces a time matrix for after infection phase (i)
    vtemp = numpy.empty([Ny,Nx,3],dtype = float)#Produces a matrix that will be filled with the amount virus above each cell
    th = numpy.empty([Ny,Nx],dtype = float) #Produces a time matrix hor healthy cells (t)
    ut = numpy.empty([Ny,Nx],dtype = float) #Produces a univeral time matrix (ut)
    tsmatrix = numpy.empty([Ny,Nx],dtype = float)

    ############################################################################
    # Filling matrices
    ############################################################################
    cells.fill("0")
    ecl.fill(0.0)
    inf.fill(0.0)
    vtemp.fill(0.0)
    th.fill(0.0)
    ut.fill(0.0)
    tsmatrix.fill(timestep) #Produces a matrix filled with value of time step

    ############################################################################
    # Deletes past versions of cells_over_time and virus_over_time
    #             and fills them with initial values
    ############################################################################
    if(BigIndex == 0):
        with open(os.path.join(Path_to_Folder,"cells_over_time.txt"),'w') as outfile:
            for i in range(2*NumberOfLayers-1):
                print(", ".join(["{}".format(LocationData[i][j]) for j in range(2*NumberOfLayers-1)]),file = outfile)
                for j in range(2*NumberOfLayers-1):
                    cells[i][j] = LocationData[i][j]
        with open(os.path.join(Path_to_Folder,"virus_over_time.txt"),'w') as outfile:
            for cell in cells:
                print("".join(["{},".format(0.0) for i in cell]),file = outfile)
    else:
        for i in range(2*NumberOfLayers-1):
            for j in range(2*NumberOfLayers-1):
                cells[i][j] = LocationData[i][j]

    ############################################################################
    #   Infects a random cell, now seen as (e)
    ############################################################################
    NumberOfInfectedCellsCount = 0
    while NumberOfInfectedCellsCount < Ni:
        randx = numpy.random.randint(0,Nx)
        randy = numpy.random.randint(0,Ny)
        if (cells[randy,randx] != "0") and (cells[randy,randx] == "h"):
            cells[randy,randx] = "e"
            inf[randy,randx] = Ti()
            NumberOfInfectedCellsCount = NumberOfInfectedCellsCount + 1
    #cells[int(NumberOfLayers-1),int(NumberOfLayers-1)] = "i" #Only the center cell
    #cells.fill("i")                                           #All the cells

    ############################################################################
    # Since the size of LocationVirus is constant, we     #
    # compute it first and refer later  #
    ############################################################################
    IndexVirus = numpy.where(cells != "0")
        #IndexVirus is a list of arrays, in this case two arrays
    LocationVirus = numpy.vstack(IndexVirus)
        #LocationVirus is a matirx
    NumberVirus = len(IndexVirus[0])
        #NumberVirus is a number

    IndexDead = numpy.where(cells == "d")
        #IndexDead is a list of arrays, in this case two arrays
    NumberDead = len(IndexDead[0])
        #NumberDead is a number

    AmountOfVirus = 0

    with open(os.path.join(Path_to_Folder,"PerTimeStep.txt"),'w') as outfile:
        print(0,",",
              int(NumberOfCells),",",
              0,",",
              0,",",
              str(NumberDead),",",
              0.0,",",
              file = outfile)

    ############################################################################
    #   Checks for Heisenberg status of viral diffusion            #
    ############################################################################
    if (D*timestep/(deltxprime**2) > 0.5) is True:
        print (D*timestep/(deltxprime**2))
        sys.exit("CHANGE PARAMETERS TO FIT DIFFUSION LIMITS. VALUE MUST BE UNDER 0.5. VALUE SHOWN ABOVE")

    ###################################################################
    #                                                                 #
    #                        Runs simulation                          #
    #                                                                 #
    ###################################################################
    NumberofTimeSteps = endtime/timestep
    NumberofSavedTimeSteps = int(NumberofTimeSteps/Save)
    adaptedtimestep = timestep #variable time step
    adaptedtimestepcount = 1
    timestepcount = 0    # equal to the number of ts elapsed
    while (timestepcount < NumberofTimeSteps-1):
        if (NumberDead == NumberOfCells):
            cell2cell = 0
            freecell = 0
        if (AmountOfVirus < 10**1) and (NumberDead == NumberOfCells):
            break
        #####################################
        #
        #               Begin
        #
        #####################################
        #####################################
        #       The Healthy Cells' time     #
        #####################################
        IndexHealthy = numpy.where(cells == 'h')
            #IndexHealthy is a list of arrays, in this case two arrays
        LocationHealthy = numpy.vstack(IndexHealthy)
            #LocationHealthy is a matrix
        NumberHealthy = len(IndexHealthy[0])
            #NumberHealthy is a number
        if NumberHealthy != 0:
            for j in range(NumberHealthy):
                [Row,Column] = LocationHealthy[:,j]
                    #Row is the row location of for a cell
                    #Column is the column location for a cell
                th[Row,Column] = th[Row,Column] + timestep
                    #"th" is the time matrix for healthy cells
                    #"ts" is the time step for the model
        #####################################
        #    Eclipse phase -> Infection     #
        #####################################
        IndexEclipse = numpy.where(cells == 'e')
            #IndexEclipse is a list of arrays, in this case two arrays
        LocationEclipse = numpy.vstack(IndexEclipse)
            #LocationEclpise is a matirx
        NumberEclipse = len(IndexEclipse[0])
            #NumberEclipse is a number
        if NumberEclipse != 0:
            for j in range(NumberEclipse):
                [Row,Column] = LocationEclipse[:,j]
                    #Row is the row location of for a cell
                    #Column is the column location for a cell
                if (ecl[Row,Column] + th[Row,Column]) < ut[Row,Column]:
                    cells[Row,Column] = 'i'
                    inf[Row,Column] =Ti()
                        #"ecl" is the time matrix for after eclipse phase
                        #"th" is the time matrix for healthy cells
                        #"ut" is the univeral time matrix
                        #"cells" is the matrix of cells
                        #"inf" is the time matrix for after infection phase
        #####################################
        #       Infection spreads &         #
        #####################################
        if cell2cell == 1:
            IndexInfected = numpy.where(cells == 'i')
                #IndexInfected is a list of arrays, in this case two arrays
            LocationInfected = numpy.vstack(IndexInfected)
                #LocationInfected is a matirx
            NumberInfected = len(IndexInfected[0])
                 #NumberInfected is a number
            if NumberInfected != 0:
                for j in range(NumberInfected):
                    [Row,Column] = LocationInfected[:,j]

                    #unless dictated otherwise, all cells around target cell exist (binary toggle)
                    AboveRowExists = 1
                    LeftColumnExists = 1
                    BelowRowExists = 1
                    RightColumnExists = 1

                    AboveRow = Row-1   #row coordinate above cell
                    LeftColumn = Column-1   #column coordinate left of cell
                    BelowRow = Row+1   #row coordinate below cell
                    RightColumn = Column+1   #column coordinate right of cell

                    if AboveRow < 0:         #if the cell one row up doesn't exist, it's taken out of the equation
                        AboveRowExists = 0
                        AboveRow = 0
                    if LeftColumn < 0:         #if the cell one column to the left doesn't exist, it's taken out of the equation
                        LeftColumnExists = 0
                        LeftColumn = 0
                    if BelowRow > Ny-1:      #if the cell one row down doesn't exist, it's taken out of the equation
                        BelowRowExists = 0
                        BelowRow = 0
                    if RightColumn > Nx-1:      ##if the cell one column to the right doesn't exist, it's taken out of the equation
                        RightColumnExists = 0
                        RightColumn = 0

                    if PU1()<probi*timestep:
                        if LeftColumnExists == 1 and cells[Row,LeftColumn] != "0":
                            if cells[Row,LeftColumn] == 'h':
                                cells[Row,LeftColumn] = 'e'
                                ecl[Row,LeftColumn] = Te()

                        if RightColumnExists == 1 and cells[Row,RightColumn] != "0":
                            if cells[Row,RightColumn] == 'h':
                                cells[Row,RightColumn] = 'e'
                                ecl[Row,RightColumn] = Te()

                        if AboveRowExists == 1 and cells[AboveRow,Column] != "0":
                            if cells[AboveRow,Column] == 'h':
                                cells[AboveRow,Column] = 'e'
                                ecl[AboveRow,Column] = Te()

                        if BelowRowExists == 1 and cells[BelowRow,Column] != "0":
                            if cells[BelowRow,Column] == 'h':
                                cells[BelowRow,Column] = 'e'
                                ecl[BelowRow,Column] = Te()

                        if AboveRowExists == 1 and RightColumnExists == 1 and cells[AboveRow,RightColumn] != "0":
                            if cells[AboveRow,RightColumn] == 'h':
                                cells[AboveRow,RightColumn] = 'e'
                                ecl[AboveRow,RightColumn] = Te()

                        if BelowRowExists == 1 and LeftColumnExists == 1 and cells[BelowRow,LeftColumn] != "0":
                            if cells[BelowRow,LeftColumn] == 'h':
                                cells[BelowRow,LeftColumn] = 'e'
                                ecl[BelowRow,LeftColumn] = Te()




        #####################################
        #       Prints status of cells      #
        #                                   #
        #####################################
        if (timestepcount % Save == 0) and (BigIndex == 0):
            with open(os.path.join(Path_to_Folder,"cells_over_time.txt"),'a') as outfile:
                for cell in cells:
                    print(", ".join(["{}".format(i) for i in cell]),file = outfile)

        #####################################
        #            Virus Spreads          #
        #####################################
        for j in range(NumberVirus):
            [Row,Column] = LocationVirus[:,j]
                #Row is the row location of for a cell
                #Column is the column location for a cell

            #unless dictated otherwise, all cells around target cell exist (binary toggle)
            AboveRowExists = 1
            LeftColumnExists = 1
            BelowRowExists = 1
            RightColumnExists = 1

            AboveRow = Row-1   #row coordinate above cell
            LeftColumn = Column-1   #column coordinate left of cell
            BelowRow = Row+1   #row coordinate below cell
            RightColumn = Column+1   #column coordinate right of cell

            if cells[Row,Column] == 'i':
                rho2 = rho
            else:
                rho2 = 0
            #where rho2 is a placeholder variable

            if AboveRow < 0:         #if the cell one row up doesn't exist, it's taken out of the equation
                AboveRowExists = 0
                AboveRow = Row
            if LeftColumn < 0:         #if the cell one column to the left doesn't exist, it's taken out of the equation
                LeftColumnExists = 0
                LeftColumn = Column
            if BelowRow > Ny-1:      #if the cell one row down doesn't exist, it's taken out of the equation
                BelowRowExists = 0
                BelowRow = Row
            if RightColumn > Nx-1:      #if the cell one column to the right doesn't exist, it's taken out of the equation
                RightColumnExists = 0
                RightColumn = Column

            if cells[AboveRow,Column] == "0":
                AboveRow = Row
            if cells[AboveRow,RightColumn] == "0":
                AboveRow = Row
                RightColumn = Column
            if cells[Row,RightColumn] == "0":
                RightColumn = Column
            if cells[BelowRow,Column] == "0":
                BelowRow = Row
            if cells[Row,LeftColumn] == "0":
                LeftColumn = Column
            if cells[BelowRow,LeftColumn] == "0":
                BelowRow = Row
                LeftColumn = Column

            NNN = (vtemp[AboveRow,Column,1] + vtemp[AboveRow,RightColumn,1] +
                       vtemp[Row,RightColumn,1] + vtemp[BelowRow,Column,1] +
                   vtemp[Row,LeftColumn,1] + vtemp[BelowRow,LeftColumn,1] )

            #ALL THIS WORK DONE SO THAT THE EQUATION IS GENERALIZED FULLY BELOW#
            VirusProduced = rho2*timestep
            VirusDecay = c*vtemp[Row,Column,1]*timestep
            VirusOut = 4*Dtsx2*vtemp[Row,Column,1]
            VirusIn = 2*Dtsx2*NNN/3

            vtemp[Row,Column,2] = vtemp[Row,Column,1] + VirusProduced - VirusOut + VirusIn - VirusDecay
            if vtemp[Row,Column,2] < 1e-10:
                vtemp[Row,Column,2] = 0.
            #####################################
            #       probability of infect       #
            #        adaptive time step         #
            #####################################
            if freecell == 1:
                pinfect = vtemp[Row,Column,2]*beta*adaptedtimestep
                while (pinfect > 1.0):
                    adaptedtimestep = adaptedtimestep/2
                    pinfect = vtemp[Row,Column,2]*beta*adaptedtimestep
                    adaptedtimestepcount = adaptedtimestepcount*2
                if pinfect <= 1.0:
                    if adaptedtimestepcount != 1:
                        pinfect = vtemp[Row,Column,2]*beta*adaptedtimestep
                    while adaptedtimestepcount != 1:
                        if PU1() < pinfect:
                            if cells[Row,Column] == 'h':
                                cells[Row,Column] = 'e'
                                ecl[Row,Column] = Te()
                        adaptedtimestepcount = adaptedtimestepcount/2
                        adaptedtimestep = adaptedtimestep*2
                        pinfect = vtemp[Row,Column,2]*beta*adaptedtimestep
                    if adaptedtimestepcount == 1:
                        vtemp[Row,Column,2] = vtemp[Row,Column,1] + VirusProduced - VirusOut + VirusIn - VirusDecay
                        if PU1() < pinfect:
                            if cells[Row,Column] == 'h':
                                cells[Row,Column] = 'e'
                                ecl[Row,Column] = Te()
            adaptedtimestepcount = 1

        vtemper = vtemp[:,:,2]
        vtemp[:,:,1] = vtemper
        AmountOfVirus = numpy.sum(vtemper)
        if (timestepcount % Save == 0) and (BigIndex == 0):
            with open(os.path.join(Path_to_Folder,"virus_over_time.txt"),'a') as outfile:
                for x in vtemper:
                    print("".join(["{},".format(i) for i in x]),file=outfile)
        ##########################################
        #               kills cells              #
        ##########################################
        IndexInfected = numpy.where(cells == 'i')
            #IndexInfected is a list of arrays, in this case two arrays
        LocationInfected = numpy.vstack(IndexInfected)
            #LocationInfected is a matirx
        NumberInfected = len(IndexInfected[0])
             #NumberInfected is a number
        if NumberInfected != 0:
            for j in range(NumberInfected):
                [Row,Column] = LocationInfected[:,j]
                    #Row is the row location of for a cell
                    #Column is the column location for a cell
                if ut[Row,Column] > (inf[Row,Column] + ecl[Row,Column] + th[Row,Column]):
                    cells[Row,Column] = 'd'
                        #"ut" is the univeral time matrix
                        #"inf" is the time matrix for after infection phase
                        #"ecl" is the time matrix for after eclipse phase
                        #"th" is the time matrix for healthy cells
                        #"cells" is the matrix of cells
        IndexDead = numpy.where(cells == 'd')
            #IndexDead is a list of arrays, in this case two arrays
        NumberDead = len(IndexDead[0])
        numberDeadBefore = NumberDead
            #NumberDead is a number



















        #############################################################
        # CELL REGENERATION
        #############################################################

        locationDead = numpy.vstack(IndexDead)

        if NumberDead != 0:

            for j in range(NumberDead):

                [row, column] = locationDead[:, j]
                # row is the row location of for a cell
                # column is the column location for a cell

                if ut[row, column] > inf[row, column] + ecl[row, column] + th[row, column] + exponentialDistro():
                    ecl[row, column] = Te()
                    th[row, column] = 0
                    inf[row, column] = Ti()

                    cells[row, column] = 'h'

                    print("A regen occured")

                    locationDead = numpy.vstack(IndexDead)
                    indexDead = numpy.where(cells == 'd') # IndexDead is a list of arrays, in this case two arrays
                    NumberDead = len(IndexDead[0]) # NumberDead is a number

                    # NumberDead -= 1

                    # numberDeadAfter = NumberDead

                    totalRegenerations += 1
                    # print("{} {} {}".format(totalRegenerations, numberDeadBefore, numberDeadAfter))


































        if((timestepcount % Save) == 0):

            # day = (timestepcount * timestep)//24
            # progressInfo = str(Tc.BLUE + "{}\t\t " + Tc.GREEN + "{}\t\t " + Tc.YELLOW + "{}\t\t " + Tc.RED + "{}\t\t" + Tc.INFO + Tc.BOLD + "{}" + Tc.END).format(day, NumberHealthy, str(NumberEclipse) + "/" + str(NumberInfected), NumberDead, totalRegenerations)
            # # clearCharacters = "\r%s" # Clears current line and overwrites with updates text
            # clearCharacters = ""
            # sys.stdout.write(clearCharacters + progressInfo)

            print(int(timestepcount*timestep),"\t",
                  str(NumberHealthy),"\t\t",
                  str(NumberEclipse),"\t\t",
                  str(NumberInfected),"\t\t",
                  str(NumberDead),"\t\t",
                  str(totalRegenerations))


            with open(os.path.join(Path_to_Folder,"PerTimeStep.txt"),'a') as outfile:
                print(int(timestepcount*timestep),",",
                      str(NumberHealthy),",",
                      str(NumberEclipse),",",
                      str(NumberInfected),",",
                      str(NumberDead),",",
                      str(AmountOfVirus),",",
                      str(totalRegenerations),",",
                      file = outfile)

        ##########################################
        #   garbages unused/unreferenced memory  #
        ##########################################
        #progress bar
        # try:
        #     # if((timestepcount%(24*int(1/timestep))) == 0):
        #     #     print("Day ",(timestepcount * timestep)//24);
        #     #     print(int(timestepcount * timestep),"\t",
        #     #           str(NumberHealthy),"\t",
        #     #           str(NumberEclipse),"\t",
        #     #           str(NumberInfected),"\t",
        #     #           str(NumberDead),"\t",
        #     #           str(AmountOfVirus),"\t",
        #     #           str(totalRegenerations))
        # except:
        #     "No Progress Bar"

        gc.collect()

        #####################################
        #       The Universal Time          #
        #       is kept here (ut)           #
        #####################################
        ut = ut + tsmatrix
        timestepcount = timestepcount+1

    ################################################################
    #   Writes a file with all of our parameters/variables         #
    ################################################################
    with open(os.path.join(Path_to_Folder,"Parameters.txt"),'a') as outfile:
        print("Nx = ",str(Nx),
              "\nNy = ",str(Ny),
              "\nNt = ",str(NumberofSavedTimeSteps),
              "\nts = ",str(timestep),
              "\nInitial End Time = ",str(endtime),
              "\nActual Hours Simulated = ",str(timestepcount*timestep),
              "\nrho = ",str(rho),
              "\nD = ",str(D),
              "\ndelta x = ",str(deltxprime),
              "\nc = ",str(c),
              "\nprobability of C2C infection:",str(probi),
              "\nCells Simulated = ",str((deltxprime*2/deltx)*Nx*Ny),
              "\nSave Path: ",Path_to_Folder,
              file = outfile)
    print("BigIndex is ", BigIndex)
print("DONE")
