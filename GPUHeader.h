#include <stdio.h>
#include "cuda.h"
#include <time.h>
#include "vector_types.h"

// #define DEBUG
#ifdef DEBUG
    // #define DEBUG_MOMT_KERNEL
    // #define DEBUG_MASS_KERNEL
    // #define DEBUG_BOUNDARY_KERNEL
    // #define DEBUG_MAXAMP_KERNEL
    #ifndef DEBUG_FUNC
        #define DEBUG_FUNC
        #define ERROR 1.0e-4
        #define CHKR 720
        #define CHKC 304
        #define CHKSI 10
        #define ID_hst(row,col) col*size_hst[0] + row
        #define ID2E_hst(row,col,dim) size_hst[2]*dim + ID_hst(row,col)
        extern float* tmpout;
    #endif
#endif

#ifndef CONSTS
#define CONSTS
    #define GX  1.0e-5
    #define EPS 1.0e-10
    #define TWLVTH 0.0833333333333333333
    #define GRAV 9.807
#endif


#ifndef CUDA_CHK
    #define CUDA_CHK
        #define cudaCHK(FUN) ({\
            if ((FUN) != cudaSuccess) {\
                printf("%s in %s at line %d\n",cudaGetErrorString(cudaGetLastError()), __FILE__,__LINE__);\
                exit(EXIT_FAILURE);\
            }\
        })
        #define cudaERROR(err) ({\
            if (err != cudaSuccess) {\
                printf("error code: %d\n", err);\
                printf("%s in %s at line %d\n",cudaGetErrorString(err), __FILE__,__LINE__);\
                exit(EXIT_FAILURE);\
            }\
        })
#endif

#ifndef CUDA_KERNEL
    #define CUDA_KERNEL
    #define ID(row,col) (col)*size_dev[0] + row
    #define ID2E(row,col,dim) size_dev[2]*(dim) + ID(row,col)
    #define F2_PITCH_ACCESS(BASS, PITCH, ROW, COL) ((float2*)((char*)(BASS) + (COL)*(PITCH)) + (ROW))
#endif


#ifndef CUDA_GLOB_VAR
    extern float *Zout_hst/*, *MNout_hst */;
    extern float /**MNdat_hst, */*Zdat_hst;
    extern float *R24_hst, *R35_hst, *H_hst;
    //float *R1_hst, *R6_hst, *R11_hst;
    extern float *R_MASS_hst;

    // extern __device__ float *R35_dev;
    // extern __device__ float *R24_dev, *H_dev;
    // extern __device__ float *Z_dat_dev, *MN_dat_dev;
    extern __device__ float/* *MN_out_dev,*/ *Z_out_dev;
    extern __constant__ __device__ uint32_t size_dev[4];
    extern float *Zmax_hst;
    extern uint32_t size_hst[4];
    extern cudaDeviceProp dev_prop;


    extern texture <float2, cudaTextureType2D, cudaReadModeElementType> MNtext;
    extern float2 *MNcontainer; //temporary space to store the float2 array
    extern size_t MNpitch;
    extern float2 *MNdat_pitchedMEM_hst;
    extern float2 *MNout_pitchedMEM_hst;
    extern __device__ float2 *MNout_pitchedMEM_dev;
    extern __constant__ __device__ size_t MNpitch_dev;

    // extern __constant__ size_t textOffset;
#endif
