# GPU-comcot

This is a program that offloads the computational component of COMCOT to Nvidia GPU. Currently, the speed up achieved by parallelized code on GTX-1060 comparing to serial one on AMD-FX8150 is nearly 200X. The code is still under development to fulfill the full functionality of the original model.

---

## **Requirements**

### **software** 

>(possible supported version/development environment):

- gfortran (4.7/4.8)

- nvcc (/9.2)

- GNU make (/3.82)

> if you wish to use our sample post-processing script. A python intepretor with the following pakages are needed:

- python(3/3.6)

- numpy(1.14)

- matplotlib(2.2)

- basemap(1.1.0)

### **hardware:**

- Nvidia-Graphic-Card with Compute Capability (3.5/6.1)

## **Installation**

### clone

```shell
git clone https://github.com/HandsomeAndy/GPU-comcot.git
```

### Setup

#### GNU solftwares

- Generally, gnu `make` and `gfortran` are pre-installed on linux machines. Try `make --version` and `gfortran --version` to check if they exist. If not, use package manager to install them.    

#### CUDA

- Check whether cuda and the card are installed and worked correctly. ```nvcc --version ``` would show the version of nvidia C compiler, and ```nvidia-smi``` list the Nvidia units ready for GPGPU, if the driver and toolkits are all set. [here](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/) is the full guide to get it work.

#### Makefile

- Modify the first few lines in `Makefile` where indicating the location of your cuda library, compute capability, and name of compilers. The example `Makefile` works on CentOS whose CUDA library is installed by `yum` equiped with 6.1 compute capability Nvidia's cards.

### Install

```shell
make
```

---

## Run a Simulation

- The input file of COMCOT is called `comcot.ctl`. Modify the file to specify a simulation time, time step, fault parameters, and grid settings. Additionally, a corresponding topographic file (available at [etopo](https://www.ngdc.noaa.gov/mgg/global/)) is required for grid construction, and the path to the file should be provided in `comcot.ctl`. For example, the fresh downloaded `comcot.ctl` takes `../etopo_halk2.xyz` as the topographic file and simulate the 2011 Tōhoku tsunami.

- To run the simulation, simply execute the program `comcot`

 ```shell
 ./comcot
 ```

## Post-Processing

- The program outputs several files at specified steps during runtime. The naming policy of those files follows the original `comcot`, while using binary as files format instead of ASCII for performance reason. The sample script `plot_dat.py` can be invoked once the simulation begins. It detects outputed files automatically and processes a series of `.png` plots with respect to the tsunami waveheight. 

```shell
python plot_dat.py .
```


## Features


## **Acknowledgements**

GPU version of COMCOT was developed by Tao, Chiu at Tsunami reseach group, IHOS, NCU and the GPU codes are protected under GPL v3.0. The goal of this work is to librate, coporate ideas with the community and accelerate the development of a high-throughput tsunami warning system at a relatively low cost. Original COMCOT version can be found at [here](http://223.4.213.26/archive/tsunami/cornell/comcot_down.htm).
