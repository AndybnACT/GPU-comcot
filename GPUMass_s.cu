#include "GPUHeader.h"
#include "GPUConfig.h"

__global__ void mass_kernel(struct GPU_Layer);



extern "C" void mass_launch_(const float* Z_f, float* Z_f_complete, const float *H_f, const int *lid){
    cudaError_t err;
    clock_t st, fi;
    struct GPU_Layer *L = ldlayer(*lid);

    st = clock();
    mass_kernel <<< L->DimGridMass, DimBlockMass >>> (*L);
    cudaDeviceSynchronize();
    err = cudaGetLastError();
    cudaERROR(err);
    fi = clock();

#ifdef DEBUG
    printf("TIME SPENT ON GPU %f\n",(float)(fi-st)/CLOCKS_PER_SEC);

#endif /* DEBUG */
}

__global__ void mass_kernel(struct GPU_Layer L){
                                /*+-->+-->+---->|
                                  +-->+-->+---->|
                                  +-->+-->+---->|
                                  +-->+-->+---->|
                                  */
                                  
    const float* __restrict__ Z = L.Zdat_hst;
    const float* __restrict__ MN = L.MNdat_hst;
    const float* __restrict__ R_MASS = L.R_MASS_hst;
    const float* __restrict__ H = L.H_hst;
    const uint32_t __restrict__ *size_dev = all_size_dev[L.lid];
    float* __restrict__ Z_out_dev = L.Zout_hst;
    
    //designed for architectures whose warpsize=32
    uint32_t row = blockIdx.x*31*(blockDim.x>>5) + 31*(threadIdx.x>>5) + threadIdx.x%32;
    uint32_t col = blockIdx.y*(size_dev[1]/gridDim.y);
    uint32_t col_end = (blockIdx.y == gridDim.y-1)? size_dev[1]-1:(blockIdx.y+1)*(size_dev[1]/gridDim.y)+1;
    float h,z;
    float m, m_suf;
    float n, n_prev;
    float ztmp;
    float r1, r11;
    float r6, r6_prev;

    n_prev = MN[ID2E(row,col,1)];
    r6_prev = R_MASS[col*4+1];

    for (uint32_t i = col+1; i < col_end; i++) {
        if (threadIdx.x%32 == 0) {
            r1  = R_MASS[i*4];
            r6  = R_MASS[i*4+1];
            r11 = R_MASS[i*4+2];
        }
        __syncwarp();
        r1 = __shfl_sync(0xFFFFFFFF,r1,0);
        r6 = __shfl_sync(0xFFFFFFFF,r6,0);
        r11 = __shfl_sync(0xFFFFFFFF,r11,0);
        m = MN[ID(row,i)];
        h =  H[ID(row,i)];
        z =  Z[ID(row,i)];
        n = MN[ID2E(row,i,1)];
        m_suf = __shfl_up_sync(0xFFFFFFFF,m,1);
        if (threadIdx.x%32 != 0 && row < size_dev[0]-1) {
            ztmp = z - r1*(m-m_suf) - r11*(n*r6-n_prev*r6_prev);
            if (ztmp + h <= EPS) ztmp = -h;
            if (h <= GX || (ztmp < EPS && -ztmp < EPS) ) ztmp = 0.0;
            Z_out_dev[ID(row,i)] = ztmp;

            r6_prev = r6;
            n_prev = n;
        }
    }
}
