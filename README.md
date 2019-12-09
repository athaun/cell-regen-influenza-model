# TCU Research Apprenticeship Program - Viral Model with regeneration

## Running python version (working)
**(works best in linux using anaconda 3.7)**

Download/clone (and extract) to any directory, which will be referenced as `base/`

CD to ```base/pythonViralTrans/```
run:
```bash
$ chmod +X multithread.sh
$ ./multithread.sh
```
to run an individual test, CD to ```base/pythonViralTrans/``` and run
```bash
$ python3 2DProbabilityViralTransmission.py ***
```
*** = the value passed through the command line as a parameter for cell regeneration

## Python code requirements

* [python 3.7](https://www.python.org/downloads/) (I recommend just downloading [anaconda 3.7](https://www.anaconda.com/distribution/))
* [numpy](https://numpy.org/) (included in Anaconda)

## Running (cuda) version (work in progress):

Download/clone (and extract) to any directory, which will be referenced as `base/`

CD to ```base/CViralTrans/```
run:
```bash
$ nvcc 1.6-CUDA-2DViralTransmission.cu -o vt.out && ./vt.out
```

## C code requirements

* [Nvidia Cuda Toolkit](https://developer.nvidia.com/cuda-downloads)
* GCC (C Compiler)
