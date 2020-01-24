#include <stdio.h>
#include "cuda.h"
#include <time.h>

//#define DEBUG
#ifdef DEBUG
    #ifndef DEBUG_FUNC
        #define DEBUG_FUNC
        #define ERROR 1.0e-4
        #define CHKR 720
        #define CHKC 304
        #define CHKSI 10
        #define ID_hst(row,col) col*size_hst[0] + row
        //#define ID2E_hst(row,col,dim) size_hst[2]*dim + ID(row,col)
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
#endif

#define MAX_LAYERS 20   
struct GPU_Layer {
    int lid;
    int plid;
    int lvl;
    float *Zdat_hst, *MNdat_hst;
    float *Zout_hst, *MNout_hst;
    float *R24_hst, *R35_hst, *H_hst;
    float *R_MASS_hst;
    float *Zmax_hst;
    uint32_t size_hst[4];
    
    dim3 DimGridMomt_MN;
    dim3 DimGridMomt;
    dim3 DimGridMass;
    dim3 DimGridOpenBD_LR;
    dim3 DimGridOpenBD_TB;
    size_t GridMaxAmp;
};
extern struct GPU_Layer Layer_struct[MAX_LAYERS];
static inline struct GPU_Layer* ldlayer(int lid){
    if (lid >= MAX_LAYERS || lid < 0) {
        printf("invalid layer id %d\n", lid);
        return NULL;
    }
    return Layer_struct + lid;
}

#ifndef CUDA_GLOB_VAR
    extern float *Zout_hst, *MNout_hst;
    extern float *MNdat_hst, *Zdat_hst;
    extern float *R24_hst, *R35_hst, *H_hst;
    //float *R1_hst, *R6_hst, *R11_hst;
    extern float *R_MASS_hst;

    // extern __device__ float *R35_dev;
    // extern __device__ float *R24_dev, *H_dev;
    // extern __device__ float *Z_dat_dev, *MN_dat_dev;
    extern __device__ float *MN_out_dev, *Z_out_dev;
    extern __constant__ __device__ uint32_t all_size_dev[MAX_LAYERS][4];
    // extern texture <float, cudaTextureType2D, cudaReadModeElementType> ZtexRef;
    // extern texture <float, cudaTextureType2D, cudaReadModeElementType> MtexRef;
    // extern texture <float, cudaTextureType2D, cudaReadModeElementType> NtexRef;
    extern float *Zmax_hst;
    extern uint32_t all_size[MAX_LAYERS][4]; // mirror of all_size_dev
    extern cudaDeviceProp dev_prop;
#endif
