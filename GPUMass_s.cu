#include "GPUHeader.h"
#include "GPUConfig.h"

extern "C" void mass_launch_(const float*, float*, const float*);
__global__ void mass_kernel(const float* __restrict__,const float* __restrict__,
                    const float* __restrict__,const float* __restrict__);
__global__ void mass_kernel_pitched(const float* __restrict__,const float2* __restrict__,
                    const float* __restrict__,const float* __restrict__);

__global__ void mass_kernel_pitched_(const float2* __restrict__, const float2* __restrict__,
                            const float* __restrict__);


extern "C" void mass_launch_(const float* Z_f, float* Z_f_complete, const float *H_f){
    cudaError_t err;
    clock_t st, fi;

    // cudaCHK( cudaMemcpy(H_hst, H_f, size_hst[3], cudaMemcpyHostToDevice) );
    //cudaCHK( cudaMemcpy(Zdat_hst, Z_f, size_hst[3], cudaMemcpyHostToDevice) );
    // printf("MASS_LAUNCH\n");
    st = clock();
    mass_kernel_pitched <<< DimGridMass, DimBlockMass >>> (Zdat_hst, MNdat_pitchedMEM_hst, R_MASS_hst, H_hst);
    cudaDeviceSynchronize();
    err = cudaGetLastError();
    cudaERROR(err);
    fi = clock();

    st = clock();
    mass_kernel_pitched_ <<< DimGridMass, DimBlockMass >>> (ZHdat_pitchedMEM_hst, MNdat_pitchedMEM_hst, R_MASS_hst);
    cudaDeviceSynchronize();
    err = cudaGetLastError();
    cudaERROR(err);
    fi = clock();




    #ifdef DEBUG_MASS_KERNEL

        st = clock();
        mass_kernel <<< DimGridMass, DimBlockMass >>> (Zdat_hst, MNdat_hst, R_MASS_hst, H_hst);
        cudaDeviceSynchronize();
        err = cudaGetLastError();
        cudaERROR(err);
        fi = clock();
        // printf("TIME SPENT ON GPU %f\n",(float)(fi-st)/CLOCKS_PER_SEC);
        // printf("printing information for debugging\n" );
        // cudaCHK( cudaMemcpy(tmpout, Zout_hst, size_hst[3], cudaMemcpyDeviceToHost) );
        // for (size_t i = 0; i < size_hst[2]; i++) {
        //     if (abs(tmpout[i] - Z_f_complete[i]) > ERROR) {
        //         printf("Z[%d] Z_cu:%e Z_f:%e %e\n", i, tmpout[i], Z_f_complete[i], tmpout[i] - Z_f_complete[i]);
        //     }
        // }
    #else
        // cudaCHK( cudaMemcpy(Z_f_complete, Zout_hst, size_hst[3], cudaMemcpyDeviceToHost) );
    #endif
}

__global__ void mass_kernel(const float* __restrict__ Z, const float* __restrict__ MN,
                            const float* __restrict__ R_MASS, const float* __restrict__ H){
                                /*+-->+-->+---->|
                                  +-->+-->+---->|
                                  +-->+-->+---->|
                                  +-->+-->+---->|
                                  */
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
            #ifdef DEBUG_MASS_KERNEL
                float z_from_pitched_kernel = Z_out_dev[ID(row,i)];
                if (fabs(z_from_pitched_kernel - ztmp) > ERROR) {
                    printf("Z/Z_[%d,%d] inconsistent: Z_%e\tZ:%e\n",row,i, z_from_pitched_kernel, ztmp);
                }
            #endif
            Z_out_dev[ID(row,i)] = ztmp;

            r6_prev = r6;
            n_prev = n;
        }
    }
}

__global__ void mass_kernel_pitched(const float* __restrict__ Z, const float2* __restrict__ MN,
                            const float* __restrict__ R_MASS, const float* __restrict__ H){
                                /*+-->+-->+---->|
                                  +-->+-->+---->|
                                  +-->+-->+---->|
                                  +-->+-->+---->|
                                  */
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

    n_prev = F2_PITCH_ACCESS(MN, MNpitch_dev, row, col)->y;//MN[ID2E(row,col,1)];
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
        m = F2_PITCH_ACCESS(MN, MNpitch_dev, row, i)->x;//MN[ID(row,i)];
        n = F2_PITCH_ACCESS(MN, MNpitch_dev, row, i)->y;//MN[ID2E(row,i,1)];
        h =  H[ID(row,i)];
        z =  Z[ID(row,i)];
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

__global__ void mass_kernel_pitched_(const float2* __restrict__ ZH, const float2* __restrict__ MN,
                            const float* __restrict__ R_MASS){
                                /*+-->+-->+---->|
                                  +-->+-->+---->|
                                  +-->+-->+---->|
                                  +-->+-->+---->|
                                  */
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

    n_prev = F2_PITCH_ACCESS(MN, MNpitch_dev, row, col)->y;//MN[ID2E(row,col,1)];
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
        m = F2_PITCH_ACCESS(MN, MNpitch_dev, row, i)->x;//MN[ID(row,i)];
        n = F2_PITCH_ACCESS(MN, MNpitch_dev, row, i)->y;//MN[ID2E(row,i,1)];
        h = F2_PITCH_ACCESS(ZH, MNpitch_dev, row, i)->y;//
        z = F2_PITCH_ACCESS(ZH, MNpitch_dev, row, i)->x;//Z[ID(row,i)];
        m_suf = __shfl_up_sync(0xFFFFFFFF,m,1);
        if (threadIdx.x%32 != 0 && row < size_dev[0]-1) {
            ztmp = z - r1*(m-m_suf) - r11*(n*r6-n_prev*r6_prev);
            if (ztmp + h <= EPS) ztmp = -h;
            if (h <= GX || (ztmp < EPS && -ztmp < EPS) ) ztmp = 0.0;
            ANSCHK(Z_out_dev[ID(row, i)], ztmp);
            F2_PITCH_ACCESS(ZHout_pitchedMEM_dev, MNpitch_dev, row, i)->x = ztmp;// Z_out_dev[ID(row,i)] = ztmp;
            // if (row == 2807 && i == 1477) {
            //     printf("Mass out=%e\narguments:%e\t%e\t%e\t%e\t%e\n%e\t%e\t%e\t%e\t%e\n----------------------\n",ztmp,r1,r6,r6_prev,r11,h,m,m_suf,n,n_prev,z );
            // }
            r6_prev = r6;
            n_prev = n;
        }
    }
}
