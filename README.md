<p align="center">
  <img src="/Splash.png" alt="Logo">

  <h2 align="center">Investigating the Effects of Cellular regeneration on viral transmission Using an Agent-Based Model</h2>

  <p align="center">
    TCU Research Apprenticeship program
    <br />
    <img src="https://img.shields.io/badge/Made%20using-NVidia%20Cuda-brightgreen">
    <img src="https://img.shields.io/badge/Made%20Using-Python%203-yellow">
    <img src="https://img.shields.io/badge/Version-2.0-blue">
    <hr>
    Treating viruses and infections is important for the general welfare. Everyone gets sick and illness negatively affects all aspects of life. Most viral infections can last for weeks, even influenza (the flu). During infections, healthy cells can grow in order to replenish the cells dying from the virus. Past viral models, especially those for influenza, tend to ignore cellular regeneration – even when it occurs in short duration infections like the flu. This research accounts for cellular regeneration, using an agent-based framework, and varies the regeneration rate in order to understand how this regeneration affects viral infections. The model used represents virus infections and spread in a two-dimensional layer of cells in order to generate total virus over time graphs for corresponding regeneration rates.
  </p>
</p>



## Running C simulation:
### Single run
Download/clone (and extract) to any directory, which will be referenced as `base/`

cd to ```base/CViralTrans/viralTransv2.0```
run:
```bash
$ nvcc ViralTransmission.cu -o program.out && ./program.out ***
```
*** = the value passed through the command line as a parameter for cell regeneration
### Multiple runs (linux/unix only)
cd to ```base/CViralTrans/viralTransv2.0```
run:
```bash
$ chmod +x runGPU.sh
$ ./runGPU.sh
```
> Warning: Depending on the run parameters, log files can be in excess of 30GB!

## Library and API requirements:

**C code requirements:**

* [Nvidia Cuda Toolkit](https://developer.nvidia.com/cuda-downloads)
* GCC (C Compiler)
* Cuda capable GPU with compute capability > 1.0

**Python graphing code requirements:**

* [python 3.7](https://www.python.org/downloads/) (I recommend just downloading [anaconda 3.7](https://www.anaconda.com/distribution/))
* [numpy](https://numpy.org/) (included in Anaconda)
* [matplotlib](https://matplotlib.org/)
