# TCU Research Apprenticeship Program - Viral Model with regeneration

## Running python simulation:
**(works best in linux using anaconda 3.7)**

Download/clone (and extract) this repository to any directory, the extracted directory will be referenced as `base/`

cd to ```base/pythonViralTrans/``` and run:
```bash
$ chmod +X multithread.sh
$ ./multithread.sh
```
to run an individual test, cd to ```base/pythonViralTrans/``` and run
```bash
$ python3 2DProbabilityViralTransmission.py ***
```
*** = the value passed through the command line as a parameter for cell regeneration

## Running C (cuda) simulation:

Download/clone (and extract) to any directory, which will be referenced as `base/`

cd to ```base/CViralTrans/```
run:
```bash
$ nvcc ViralTransmission v1.6.cu -o program.out && ./program.out
```

## Creating graphs:
### To create graphs of the simulation over time, including healthy cells vs time, virus vs time, dead cells vs time, and more, follow the instructions below.
**1) Find the path to ViralModel/:**

* In your file manager, copy the path to `ViralModel/`, for example `~/Documents/pythonViralTrans/ViralModel/`.

* Open `createGraphs.py` in any text editor, and paste the path into the variable named `OuterPath`.

You will need to adjust the parameters in `createGraphs.py` to match the parameters set in the simulation code, you can follow either of the below instruction sets to go about different ways of finding these values, but I highly recommend you start with the first set especially if you ran the python simulation as it is easier to follow, and all you need; however if you ran the C simulation you will still need to follow some of the steps outlined in #3.

**2) Copy values from simulation code:**<br>
<br>If you ran the python simulation, all of the values you need can be found near the top of the simulation programs in a section labeled "simulation parameters" (*lines 73-81 in the python code, and lines 41-48 in the C code*). If you ran the C code most of the values can be found here, except for MOI which must be determined by using step 3.2.
<br><br>

**3) Find parameters for createGraphs.py based on the folders in ViralModel/ (Hard way):** <br><br>
An example of a folder name could look like this: `ViralModel/607_0-FREECELL_-2.0-MOI_0.01-RP/`
* `NumberOfLayers` can be determined by the first value, `607`.

* The next parameter, `cellToCell`, can be determined by the string following `607_0-`; if that string says "FREECELL", set `cellToCell` to `False`, if it says "CELL2CELL", set it to `True`.

* MOI is the next value after `FREECELL_`, in the case of our example `-2.0`, but it may be different, and their may be multiple folders with different values of MOI, so add them all into the array/list

* `regenParams` or the regeneration parameters may be different on each folder, and they can be found between `MOI_` and `-RP`. add each value to the array/list, making sure that they are floating point values.

* `NumberOfRuns` must be determined by checking the simulation code, in the python version is can be found on line `79` of `2DProbabilityViralTransmission.py`, and line `48` of `ViralTransmission v1.6.cu`.

## Library and API requirements:

**C code requirements:**

* [Nvidia Cuda Toolkit](https://developer.nvidia.com/cuda-downloads)
* GCC (C Compiler)

**Python code requirements (including graphing code):**

* [python 3.7](https://www.python.org/downloads/) (I recommend just downloading [anaconda 3.7](https://www.anaconda.com/distribution/))
* [numpy](https://numpy.org/) (included in Anaconda)
* [matplotlib](https://matplotlib.org/) (used for the graphing code, not required to run simulation)
