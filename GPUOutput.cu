#include "GPUHeader.h"
#include "GPUConfig.h"

extern "C" void maxamp_launch_(void);
__global__ void maximum_recorder_kernel(float *, const float * __restrict__ , const uint32_t);
extern "C" void cuda_getz_(float *);
extern "C" void cuda_getmn_(float *,float *);
extern "C" void cuda_getzmax_(float *);


extern "C" void maxamp_launch_(void){
    clock_t st, fi;
    cudaError_t err;

    st = clock();
    maximum_recorder_kernel <<< GridMaxAmp, MAXAMP_BLOCK >>> (Zmax_hst, Zout_hst, size_hst[2]);
    cudaDeviceSynchronize();
    err = cudaGetLastError();
    cudaERROR(err);
    fi = clock();

    #ifdef DEBUG_MAXAMP_KERNEL
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


extern "C" void cuda_getz_(float *Z_f) {
    cudaCHK( cudaMemcpy(Z_f, Zout_hst, size_hst[3], cudaMemcpyDeviceToHost) );
}

extern "C" void cuda_getmn_(float *M_f, float *N_f) {
    cudaCHK( cudaMemcpy2D(MNcontainer, size_hst[0]*sizeof(float2), MNout_pitchedMEM_hst, MNpitch, size_hst[0]*sizeof(float2), size_hst[1], cudaMemcpyDeviceToHost) );
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
