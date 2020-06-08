<p align="center">
  <img src="/Splash.png" alt="Logo">

  <h2 align="center">Investigating the Effects of Cellular regeneration on viral transmission Using an Agent-Based Model</h2>

  <p align="center">
    TCU Research Apprenticeship program
    <br />
    
  </p>
</p>



# Investigating the Effects of Cellular regeneration on viral transmission Using an Agent-Based Model
### TCU Research Apprenticeship Program

Treating viruses and infections is important for the general welfare. Everyone gets sick and illness negatively affects all aspects of life. Most viral infections can last for weeks, even influenza (the flu). During infections, healthy cells can grow in order to replenish the cells dying from the virus. Past viral models, especially those for influenza, tend to ignore cellular regeneration – even when it occurs in short duration infections like the flu. This research accounts for cellular regeneration, using an agent-based framework, and varies the regeneration rate in order to understand how this regeneration affects viral infections. The model used represents virus infections and spread in a two-dimensional layer of cells in order to generate total virus over time graphs for corresponding regeneration rates.

## Running C simulation:
Note: The cuda version of this simulation works, but does not include cellular regeneration yet.

Download/clone (and extract) to any directory, which will be referenced as `base/`

cd to ```base/CViralTrans/viralTransv1.8```
run:
```bash
$ nvcc ViralTransmission.cu -o program.out && ./program.out ***
```
*** = the value passed through the command line as a parameter for cell regeneration

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

## Library and API requirements:

**C code requirements:**

* [Nvidia Cuda Toolkit](https://developer.nvidia.com/cuda-downloads)
* GCC (C Compiler)

**Python code requirements (including graphing code):**

* [python 3.7](https://www.python.org/downloads/) (I recommend just downloading [anaconda 3.7](https://www.anaconda.com/distribution/))
* [numpy](https://numpy.org/) (included in Anaconda)
* [matplotlib](https://matplotlib.org/) (used for the graphing code, not required to run simulation)
