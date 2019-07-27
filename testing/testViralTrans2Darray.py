environment = {
    "infectionRate": 1,  # 1 cell per cycle
    "cycleRate": 1,      # 1 cycle per day
    "simulatedDays": 7,
}
cellsPerDay = environment["infectionRate"] * environment["simulatedDays"]  # cells infected per day

import math

numCells = 0

class Cell:

    def __init__(self, id):
        self.infected = False
        self.dead = False
        self.hasNeighbors = False
        self.id = id

    def printId(self):
        print(self.id)

    def infect(self):
        self.infected = True


cells = []
def newCell(id):
    cells.append(Cell(id))


def createTestingEnv(numberOfCells):
    # numberOfCells = numberOfCells
    global numCells
    numCells = numberOfCells
    # creates a simulated cell tissue
    global env
    env = []
    for x in range(numberOfCells):
        env.append([])
        for y in range(numberOfCells):
            env[x].append(0)
            newCell(x * numberOfCells + y)


createTestingEnv(10 ** 5)  # creates a 2d array with 1 million indices

for a in cells:
    a.printId()

# print()
