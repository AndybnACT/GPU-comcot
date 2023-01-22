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
    maximum_recorder_kernel <<< L->GridMaxAmp, MAXAMP_BLOCK >>> (L->Zmax_hst, L->Zout_hst, L->l_size[2]);
    cudaDeviceSynchronize();
    err = cudaGetLastError();
    cudaERROR(err);
    fi = clock();

    #ifdef DEBUG
        printf("GPU MAX_AMP\n" );
        printf("TIME SPENT ON GPU %f\n",(float)(fi-st)/CLOCKS_PER_SEC);
    #endif
}

extern "C" void gcomcot_get_z_(int *lid, float *l1, float *l2, 
                                         float *l3, float *l4,
                                         int *i, int *j)
{
    struct GPU_Layer *L = ldlayer(*lid);
    float tmp[2];
    cudaCHK( cudaMemcpy(tmp, L->Zout_hst + L->l_size[0] * (*j - 1) + (*i - 1),
                         2 * sizeof(float), cudaMemcpyDeviceToHost) );
    *l1 = tmp[0];
    *l2 = tmp[1];
    cudaCHK( cudaMemcpy(tmp, L->Zout_hst + L->l_size[0] * (*j) + (*i - 1),
                         2 * sizeof(float), cudaMemcpyDeviceToHost) );
    *l3 = tmp[0];
    *l4 = tmp[1];
    return;
}

__global__ void maximum_recorder_kernel(float *max_array, const float* __restrict__ array, const uint32_t array_size){
    uint32_t id=blockIdx.x*blockDim.x + threadIdx.x;
    if (id < array_size) {
        float current_max_value = max_array[id];
        if (current_max_value >= array[id]) return;
        else max_array[id] = array[id];
    }
}

extern "C" void cuda_getz1_(float *Z_f, int *lid) {
    struct GPU_Layer *L = ldlayer(*lid);
    cudaCHK( cudaMemcpy(Z_f, L->Zdat_hst, L->l_size[3], cudaMemcpyDeviceToHost) );
}

extern "C" void cuda_setz1_(float *Z_f, int *lid) {
    struct GPU_Layer *L = ldlayer(*lid);
    cudaCHK( cudaMemcpy(L->Zdat_hst, Z_f, L->l_size[3], cudaMemcpyHostToDevice) );
}

extern "C" void cuda_getz_(float *Z_f, int *lid) {
    struct GPU_Layer *L = ldlayer(*lid);
    cudaCHK( cudaMemcpy(Z_f, L->Zout_hst, L->l_size[3], cudaMemcpyDeviceToHost) );
}

extern "C" void cuda_setz_(float *Z_f, int *lid) {
    struct GPU_Layer *L = ldlayer(*lid);
    cudaCHK( cudaMemcpy(L->Zout_hst, Z_f, L->l_size[3], cudaMemcpyHostToDevice) );
}

extern "C" void cuda_getmn1_(float *M_f, float *N_f, int *lid) {
    struct GPU_Layer *L = ldlayer(*lid);
    cudaCHK( cudaMemcpy(M_f, L->MNdat_hst, L->l_size[3], cudaMemcpyDeviceToHost) );
    cudaCHK( cudaMemcpy(N_f, L->MNdat_hst+L->l_size[2], L->l_size[3], cudaMemcpyDeviceToHost) );
}

extern "C" void cuda_setmn1_(float *M_f, float *N_f, int *lid) {
    struct GPU_Layer *L = ldlayer(*lid);
    cudaCHK( cudaMemcpy(L->MNdat_hst, M_f, L->l_size[3], cudaMemcpyHostToDevice) );
    cudaCHK( cudaMemcpy(L->MNdat_hst+L->l_size[2], N_f, L->l_size[3], cudaMemcpyHostToDevice) );
}

extern "C" void cuda_getmn_(float *M_f, float *N_f, int *lid) {
    struct GPU_Layer *L = ldlayer(*lid);
    cudaCHK( cudaMemcpy(M_f, L->MNout_hst, L->l_size[3], cudaMemcpyDeviceToHost) );
    cudaCHK( cudaMemcpy(N_f, L->MNout_hst+L->l_size[2], L->l_size[3], cudaMemcpyDeviceToHost) );
}

extern "C" void cuda_setmn_(float *M_f, float *N_f, int *lid) {
    struct GPU_Layer *L = ldlayer(*lid);
    cudaCHK( cudaMemcpy(L->MNout_hst, M_f, L->l_size[3], cudaMemcpyHostToDevice) );
    cudaCHK( cudaMemcpy(L->MNout_hst+L->l_size[2], N_f, L->l_size[3], cudaMemcpyHostToDevice) );
}

extern "C" void cuda_getzmax_(float *Zmax_f, int *lid) {
    struct GPU_Layer *L = ldlayer(*lid);
    cudaCHK( cudaMemcpy(Zmax_f, L->Zmax_hst, L->l_size[3], cudaMemcpyDeviceToHost) );
}

// void cuda_getmnmax_(float *Mmax_f, float *Nmax_f) {
//     cudaCHK( cudaMemcpy(Mmax_f, MNmax_hst, l_size[3], cudaMemcpyDeviceToHost) );
//     cudaCHK( cudaMemcpy(Nmax_f, MNmax_hst+l_size[2], l_size[3], cudaMemcpyDeviceToHost) );
// }
