#include <stdio.h>
#include "debug_option.h"
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
    #define ID_MATRIX(row,col,row_size) ((col)*(row_size) + (row)) 
#endif

#define MAX_LAYERS 20   
struct GPU_Layer {
    int lid;
    int plid;
    int lvl;
    int rel_size;
    int rel_time;
    int corner[4];
    float *Zdat_hst, *MNdat_hst;
    float *Zout_hst, *MNout_hst;
    float *R24_hst, *R35_hst, *H_hst;
    float *R_MASS_hst;
    float *Zmax_hst;
    float *yflux;
    float *xflux;
    uint32_t l_size[4];
    float grx;
    float gry;
    
    dim3 DimGridJnz;
    dim3 DimBlockJnz;
    dim3 DimGridMomt_MN;
    dim3 DimGridMomt;
    dim3 DimGridMass;
    dim3 DimGridOpenBD_LR;
    dim3 DimGridOpenBD_TB;
    dim3 DimGrid_JNQV;
    dim3 DimGrid_JNQH;
    size_t GridMaxAmp;
    struct GPU_Layer *child;
    struct GPU_Layer *sibling;
};
extern struct GPU_Layer Layer_struct[MAX_LAYERS];

extern struct GPU_Layer *Root_Grid;

static inline struct GPU_Layer* ldlayer(int lid){
    if (lid >= MAX_LAYERS || lid < 0) {
        printf("invalid layer id %d\n", lid);
        return NULL;
    }
    return Layer_struct + lid;
}

extern "C" void mass_launch_(const float* Z_f, float* Z_f_complete, const float *H_f, const int *lid);
extern "C" void momt_launch_(float*,float*,float*, const int*);
extern "C" void cuda_update_layer_(int *lid);

#ifndef CUDA_GLOB_VAR
    extern __constant__ __device__ uint32_t all_size_dev[MAX_LAYERS][4];
    extern uint32_t all_size[MAX_LAYERS][4]; // mirror of all_size_dev
    extern cudaDeviceProp dev_prop;
#endif

#ifdef RELDIFF
static inline bool assert_diff(float base, float eval) {
    if (base == 0.0f) {
        if (eval != 0.0f)
            return true;
        else
            return false;
    }

    if (fabs((eval - base) / base) >= DIFFRATE)
        return true;
    else
        return false;
}
#else
static inline bool assert_diff(float base, float eval) {
    if (fabs(eval - base) >= TOLERANCE)
        return true;
    else
        return false;
}
#endif /* RELDIFF */

#define CMP_VAR(cmpp, varp, row, col, str,...) { \
    for (size_t j = 0; j < col; j++) {                          \
        for (size_t i = 0; i < row; i++) {                      \
            size_t id = j*(row) + i;                            \
            double diff = fabs((cmpp)[id] - (varp)[id]);        \
            if (assert_diff((cmpp)[id], (varp)[id])) {          \
                printf(str, __VA_ARGS__);                       \
                printf("(i,j)=(%zu,%zu), %f %f, diff=%lf  ",    \
                    i, j, (cmpp)[id], (varp)[id], diff);        \
                printf("%x, %x\n", *(unsigned int*)((cmpp) + id),\
                                *(unsigned int*)((varp) + id)); \
            }                                                   \
        }                                                       \
    }                                                           \
}

static inline bool ispow2(unsigned int val) {
    // returns true if val is power of 2
    return !(val & (val - 1));
} 

