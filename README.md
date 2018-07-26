# GPU-comcot

This is a program that offloads the computational component of COMCOT to Nvidia GPU. Currently, the speed up achieved by parallelized code on GTX-1060 comparing to serial one on AMD-FX8150 is nearly 200X. The code is still under development to fulfill the full functionality of the original model.

---

## **Requirements**

### **software** 

>(possible supported version/development environment):

- gfortran (4.7/4.8)

- nvcc (/9.2)

- GNU make (/3.82)

### **hardware:**

- Nvidia-Graphic-Card with Compute Capability (3.5/6.1)

---

## **Installation**

### clone

```shell
git clone https://github.com/HandsomeAndy/GPU-comcot.git
```

### Setup

#### GNU solftwares

- Generally, gnu `make` and `gfortran` are pre-installed on linux machines. Try `make --version` and `gfortran --version` to check if they exist. If not, use package manager to install them.    

#### CUDA

- Check whether cuda are installed and worked correctly. ```nvcc --version ``` would show the version of nvidia C compiler, and ```nvidia-smi``` list the Nvidia units ready for GPGPU, if the driver and toolkits are all set. [here](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/) is the full guide to get it work.

#### Makefile

- Modify the first few lines in `Makefile` where indicating the location of your cuda library, compute capability, and name of compilers. The example `Makefile` works on CentOS whose CUDA library is installed by `yum` equiped with 6.1 compute capability Nvidia's cards.

### Install

```shell
make
```

---

## Run the Simulation

- The input file of COMCOT called `comcot.ctl`. 

## Fetures



## **Acknowledgements**

GPU version of COMCOT was developed by Tao, Chiu at Tsunami reseach group, IHOS, NCU and protected under GPL v3.0. The goal of this work is to librate, coporate ideas with the community and accelerate the development of a high-throughput tsunami warning system at relatively low cost. Original COMCOT version can be found at [here](http://223.4.213.26/archive/tsunami/cornell/comcot_down.htm).
