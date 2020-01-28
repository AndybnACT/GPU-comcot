#include "GPUHeader.h"
#include "GPUConfig.h"

__global__ void momts_kernel(const float* __restrict__,const float* __restrict__,
                    const float* __restrict__,const float* __restrict__,const float* __restrict__);
__global__ void momt_kernelM(struct GPU_Layer L);
__global__ void momt_kernelN(struct GPU_Layer L);

extern "C" void momt_launch_(float *M_f, float *N_f, float *Z_f, const int *lid) {
    float* tmp;
    clock_t st, fi;
    cudaError_t err;
    struct GPU_Layer *L = ldlayer(*lid);

#ifdef DEBUG
    printf("Z_cu vs Z_f\n");
    cudaCHK( cudaMemcpy(tmpout, Zout_hst, L->l_size[3], cudaMemcpyDeviceToHost) );
    for (size_t i = 0; i < l_size[2]; i++) {
        if (abs(tmpout[i] - Z_f[i]) > ERROR) {
            printf("err\n");
        }
    }
    
#endif /* DEBUG */

    st = clock();
    momt_kernelM <<<  L->DimGridMomt_MN, DimBlockMomt_MN, 0, EXECstream[0] >>> (*L);
    momt_kernelN <<<  L->DimGridMomt_MN, DimBlockMomt_MN, 0, EXECstream[1] >>> (*L);
    cudaDeviceSynchronize();
    err = cudaGetLastError();
    cudaERROR(err);
    fi = clock();

#ifdef DEBUG
    printf("TIME SPENT ON GPU %f\n",(float)(fi-st)/CLOCKS_PER_SEC);
    
#endif /* DEBUG */

}

__global__ void momt_kernelM(struct GPU_Layer L) {
    
     const float* __restrict__ Z = L.Zout_hst;
     const float* __restrict__ MN = L.MNdat_hst;
     const float* __restrict__ R2 = L.R24_hst;
     const float* __restrict__ R3 = L.R35_hst;
     const float* __restrict__ H = L.H_hst;
     float* __restrict__ MN_out_dev = L.MNout_hst;
     const uint32_t __restrict__ *size_dev = all_size_dev[L.lid];
     
     
     uint32_t row = blockIdx.x*31*(blockDim.x>>5) + 31*(threadIdx.x>>5) + threadIdx.x%32;
     uint32_t col = blockIdx.y*(size_dev[1]/gridDim.y);
     uint32_t col_end = (blockIdx.y == gridDim.y-1)? size_dev[1]:(blockIdx.y+1)*(size_dev[1]/gridDim.y)+1;

     // memory operands ====================================
     float n_prev,        n; // prev-->stored for next loop
     float n_ip1jm1_prev, n_ip1;// n shuffle
     //Z
     float z;     //independant load
     float z_ip1; //Z shuffle
     //H
     float h;     //independant load
     float h_ip1; //H shuffle

     float m, r2; //independant loads
     float r3;    //boardcast
     //=====================================================
     // computation operands ===============================
     float tot_n;
     float xm;
     //=====================================================

     if (row < size_dev[0]) { //lower bound of threads
         n_prev = MN[ID2E(row,col,1)];
         n_ip1jm1_prev = __shfl_down_sync(0xFFFFFFFF, n_prev,1);
         if (col != 0) col+= 1; // not start at left boundary
         for (uint32_t j = col; j < col_end; j++) { //grid striding
             if (threadIdx.x%32 == 31) {
                 r3 = R3[j];
             }
             __syncwarp();
             r3    = __shfl_sync(0xFFFFFFFF,r3,31);
             n     = MN[ID2E(row,j,1)];
             h     = H[ID(row,j)];
             z     = Z[ID(row,j)];

             n_ip1 = __shfl_down_sync(0xFFFFFFFF, n,1);
             z_ip1 = __shfl_down_sync(0xFFFFFFFF, z,1);
             h_ip1 = __shfl_down_sync(0xFFFFFFFF, h,1);

             if (threadIdx.x%32 != 31 && row < size_dev[0]-1) { //lower bound of lanes & computation
                 if (h > GX && h_ip1 > GX) {
                     m     = MN[ID(row,j)];
                     r2    = R2[ID(row,j)];

                     tot_n = n + n_ip1 + n_prev + n_ip1jm1_prev;
                     xm    = m - r2*(z_ip1 - z) + r3*tot_n;
                     if (xm < EPS && xm > -EPS) xm = 0.0f;
                     MN_out_dev[ID(row,j)] = xm;

                 }else{
                     MN_out_dev[ID(row,j)] = 0.0f;
                 }
                 // if (row == 447 && j == 782) {
                 //     printf("[%d,%d]===============================[row,j][%d,%d]\n",threadIdx.x,blockIdx.x,row,j );
                 //     printf("%e\t%e\t%e\t%e\t%e\t%e\n",r2,r3,n,n_ip1,n_prev,n_ip1jm1_prev);
                 //     printf("%e\t%e\t%e\n",m,z_ip1,z);
                 // }

             }
             n_prev = n;
             n_ip1jm1_prev = n_ip1;
         }
     }
}



__global__ void momt_kernelN(struct GPU_Layer L) {
                                 
    const float* __restrict__ Z = L.Zout_hst;
    const float* __restrict__ MN = L.MNdat_hst;
    const float* __restrict__ R4 = L.R24_hst + L.l_size[2];
    const float* __restrict__ R5 = L.R35_hst + L.l_size[1];
    const float* __restrict__ H = L.H_hst;
    float* __restrict__ MN_out_dev = L.MNout_hst;
    const uint32_t __restrict__ *size_dev = all_size_dev[L.lid];
    

    uint32_t row = blockIdx.x*31*(blockDim.x>>5) + 31*(threadIdx.x>>5) + threadIdx.x%32;
    uint32_t col = blockIdx.y*(size_dev[1]/gridDim.y);
    uint32_t col_end = (blockIdx.y == gridDim.y-1)? size_dev[1]-1:(blockIdx.y+1)*(size_dev[1]/gridDim.y)+1;

    // memory operands ==================================== 'prev' variable should be kept in registers
    float m_im1, m_im1jp1_prev; // m shuffle & prev-->stored for next loop
    float m,        m_jp1_prev; // independant loas + prevload
    //Z
    float z, z_jp1_prev;       //independant load + prevload
    //H
    float h, h_jp1_prev;       //independant load + prevload

    float n, r4; //independant loads
    float r5;    //boardcast
    //=====================================================
    // computation operands ===============================
    float tot_m;
    float xn;
    //=====================================================



    if (row < size_dev[0]) { //lower bound of threads;
        m_jp1_prev = MN[ID(row,col_end)];
        m_im1jp1_prev = __shfl_up_sync(0xFFFFFFFF, m_jp1_prev,1);

        z_jp1_prev = Z[ID(row,col_end)];
        h_jp1_prev = H[ID(row,col_end)];
        // !! bug fixed: note that comparison of unsigned expression >= 0 is always true
        for (int j = (int)col_end-1; j >= (int)col; j--) { //grid striding
            if (threadIdx.x%32 == 0) {
                r5 = R5[j];
            }
            __syncwarp();
            r5    = __shfl_sync(0xFFFFFFFF,r5,0);
            m     = MN[ID(row,j)];
            h     =  H[ID(row,j)];
            m_im1 = __shfl_up_sync(0xFFFFFFFF,m,1);

            if (threadIdx.x%32 != 0 || row == 0) { //upper bound of lanes
                r4    = R4[ID(row,j)];
                z     =  Z[ID(row,j)];
                n     = MN[ID2E(row,j,1)];
                if (h > GX && h_jp1_prev > GX) {
                    if (row != 0) tot_m =  m_im1 + m_im1jp1_prev + m + m_jp1_prev;
                    else          tot_m =  2.0*(m + m_jp1_prev); // upper bound of the computation
                    xn    = n - r4*(z_jp1_prev -z) - r5*tot_m;

                    // if (row == 449 && j == 782) {
                    //     printf("[%d,%d]===============================[row,j][%d,%d]\n",threadIdx.x,blockIdx.x,row,j );
                    //     printf("%e\t%e\t%e\t%e\t%e\t%e\n",r4, r5, m, m_jp1_prev, m_im1jp1_prev, m_im1);
                    //     printf("%e\t%e\t%e\n",n,z_jp1_prev,z);
                    // }

                    if (xn < EPS && xn > -EPS) xn = 0.0f;
                    MN_out_dev[ID2E(row,j,1)] = xn;
                }else{
                    MN_out_dev[ID2E(row,j,1)] = 0.0f;
                }
                m_jp1_prev    = m;
                m_im1jp1_prev = m_im1;
                z_jp1_prev    = z;
                h_jp1_prev    = h;
            }
        }
    }



}
