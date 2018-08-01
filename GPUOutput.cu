#include "GPUHeader.h"
#include "GPUConfig.h"

extern "C" void maxamp_launch_(void);
// __global__ void maximum_recorder_kernel(float *, const float * __restrict__ , const uint32_t);
__global__ void maximum_recorder_kernel_tex(float *);
extern "C" void cuda_getz_(float *);
extern "C" void cuda_getmn_(float *,float *);
extern "C" void cuda_getzmax_(float *);


extern "C" void maxamp_launch_(void){
    clock_t st, fi;
    cudaError_t err;

    // st = clock();
    // maximum_recorder_kernel <<< GridMaxAmp, MAXAMP_BLOCK >>> (Zmax_hst, Zout_hst, size_hst[2]);
    // cudaDeviceSynchronize();
    // err = cudaGetLastError();
    // cudaERROR(err);
    // fi = clock();

    st = clock();
    maximum_recorder_kernel_tex <<< DimGridMaxAmp, DimBlockMaxAmp >>> (Zmax_hst);
    cudaDeviceSynchronize();
    err = cudaGetLastError();
    cudaERROR(err);
    fi = clock();



    #ifdef DEBUG_MAXAMP_KERNEL
        printf("GPU MAX_AMP\n" );
        printf("TIME SPENT ON GPU %f\n",(float)(fi-st)/CLOCKS_PER_SEC);
    #endif
}



// __global__ void maximum_recorder_kernel(float *max_array, const float* __restrict__ array, const uint32_t array_size){
//     uint32_t id=blockIdx.x*blockDim.x + threadIdx.x;
//     if (id < array_size) {
//         float current_max_value = max_array[id];
//         if (current_max_value >= array[id]) return;
//         else max_array[id] = array[id];
//     }
// }

__global__ void maximum_recorder_kernel_tex(float *max_array){
    uint32_t x=blockIdx.x*blockDim.x + threadIdx.x;
    uint32_t y=blockIdx.y*blockDim.y + threadIdx.y;
    float eta_now;
    float eta_max;
    if (x < size_dev[0] && y < size_dev[1]) {
        // if (tex2D<float2>(ZHtext, x, y).y > 0.0f) {
            eta_now = tex2D<float2>(ZHtext, x, y).x;
            eta_max = max_array[ID(x,y)];
            if (eta_max < eta_now) {
                max_array[ID(x,y)] = eta_now;
            }
        // }
    }
}

extern "C" void cuda_getz_(float *Z_f) {
    // cudaCHK( cudaMemcpy(Z_f, Zout_hst, size_hst[3], cudaMemcpyDeviceToHost) );
    cudaCHK( cudaMemcpy2D(ZHcontainer, size_hst[0]*sizeof(float2), ZHout_pitchedMEM_hst, LayerPitch, size_hst[0]*sizeof(float2), size_hst[1], cudaMemcpyDeviceToHost) );
    for (size_t i = 0; i < size_hst[1]; i++) {
        for (size_t j = 0; j < size_hst[0]; j++) {
            Z_f[i*size_hst[0] + j] = ZHcontainer[i*size_hst[0] + j].x;
        }
    }
}

extern "C" void cuda_getmn_(float *M_f, float *N_f) {
    cudaCHK( cudaMemcpy2D(MNcontainer, size_hst[0]*sizeof(float2), MNout_pitchedMEM_hst, LayerPitch, size_hst[0]*sizeof(float2), size_hst[1], cudaMemcpyDeviceToHost) );
    for (size_t j = 0; j < size_hst[1]; j++) {
        for (size_t i = 0; i < size_hst[0]; i++) {
            M_f[j*size_hst[0]+i] = MNcontainer[j*size_hst[0]+i].x;
            N_f[j*size_hst[0]+i] = MNcontainer[j*size_hst[0]+i].y;
        }
    }
    // cudaCHK( cudaMemcpy(M_f, MNout_hst, size_hst[3], cudaMemcpyDeviceToHost) );
    // cudaCHK( cudaMemcpy(N_f, MNout_hst+size_hst[2], size_hst[3], cudaMemcpyDeviceToHost) );
}

extern "C" void cuda_getzmax_(float *Zmax_f) {
    cudaCHK( cudaMemcpy(Zmax_f, Zmax_hst, size_hst[3], cudaMemcpyDeviceToHost) );
}

// void cuda_getmnmax_(float *Mmax_f, float *Nmax_f) {
//     cudaCHK( cudaMemcpy(Mmax_f, MNmax_hst, size_hst[3], cudaMemcpyDeviceToHost) );
//     cudaCHK( cudaMemcpy(Nmax_f, MNmax_hst+size_hst[2], size_hst[3], cudaMemcpyDeviceToHost) );
// }
