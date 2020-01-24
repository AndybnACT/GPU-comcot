#include "GPUHeader.h"
#include "GPUConfig.h"

extern "C" void maxamp_launch_(const int *);
__global__ void maximum_recorder_kernel(float *, const float * __restrict__ , const uint32_t);
extern "C" void cuda_getz_(float *, int *);
extern "C" void cuda_getmn_(float *,float *, int *);
extern "C" void cuda_getzmax_(float *, int *);


extern "C" void maxamp_launch_(const int *lid){
    clock_t st, fi;
    cudaError_t err;
    struct GPU_Layer *L = ldlayer(*lid);

    st = clock();
    maximum_recorder_kernel <<< L->GridMaxAmp, MAXAMP_BLOCK >>> (L->Zmax_hst, L->Zout_hst, L->size_hst[2]);
    cudaDeviceSynchronize();
    err = cudaGetLastError();
    cudaERROR(err);
    fi = clock();

    #ifdef DEBUG
        printf("GPU MAX_AMP\n" );
        printf("TIME SPENT ON GPU %f\n",(float)(fi-st)/CLOCKS_PER_SEC);
    #endif
}



__global__ void maximum_recorder_kernel(float *max_array, const float* __restrict__ array, const uint32_t array_size){
    uint32_t id=blockIdx.x*blockDim.x + threadIdx.x;
    if (id < array_size) {
        float current_max_value = max_array[id];
        if (current_max_value >= array[id]) return;
        else max_array[id] = array[id];
    }
}


extern "C" void cuda_getz_(float *Z_f, int *lid) {
    struct GPU_Layer *L = ldlayer(*lid);
    cudaCHK( cudaMemcpy(Z_f, L->Zout_hst, L->size_hst[3], cudaMemcpyDeviceToHost) );
}

extern "C" void cuda_getmn_(float *M_f, float *N_f, int *lid) {
    struct GPU_Layer *L = ldlayer(*lid);
    cudaCHK( cudaMemcpy(M_f, L->MNout_hst, L->size_hst[3], cudaMemcpyDeviceToHost) );
    cudaCHK( cudaMemcpy(N_f, L->MNout_hst+L->size_hst[2], L->size_hst[3], cudaMemcpyDeviceToHost) );
}

extern "C" void cuda_getzmax_(float *Zmax_f, int *lid) {
    struct GPU_Layer *L = ldlayer(*lid);
    cudaCHK( cudaMemcpy(Zmax_f, L->Zmax_hst, L->size_hst[3], cudaMemcpyDeviceToHost) );
}

// void cuda_getmnmax_(float *Mmax_f, float *Nmax_f) {
//     cudaCHK( cudaMemcpy(Mmax_f, MNmax_hst, size_hst[3], cudaMemcpyDeviceToHost) );
//     cudaCHK( cudaMemcpy(Nmax_f, MNmax_hst+size_hst[2], size_hst[3], cudaMemcpyDeviceToHost) );
// }
