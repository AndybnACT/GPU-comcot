#include "GPUHeader.h"
#include "GPUConfig.h"

extern "C" void momt_launch_(float*,float*,float*, int*);
__global__ void momts_kernel(const float* __restrict__,const float* __restrict__,
                    const float* __restrict__,const float* __restrict__,const float* __restrict__);
__global__ void momt_kernelM(struct GPU_Layer L);
__global__ void momt_kernelN(struct GPU_Layer L);

extern "C" void momt_launch_(float *M_f, float *N_f, float *Z_f, int *lid) {
    float* tmp;
    clock_t st, fi;
    cudaError_t err;
    struct GPU_Layer *L = ldlayer(*lid);

    #ifdef DEBUG
        printf("Z_cu vs Z_f\n");
        cudaCHK( cudaMemcpy(tmpout, Zout_hst, size_hst[3], cudaMemcpyDeviceToHost) );
        for (size_t i = 0; i < size_hst[2]; i++) {
            if (abs(tmpout[i] - Z_f[i]) > ERROR) {
                printf("Z[%d,%d] Z_cu:%e Z_f:%e %e\n", i%size_hst[0], i/size_hst[0], tmpout[i], Z_f[i], tmpout[i] - Z_f[i]);
            }
        }

    #endif

    // // kernel launch
    // st = clock();
    // momts_kernel <<< DimGridMomt, DimBlockMomt >>> (Zout_hst, MNdat_hst, R24_hst, R35_hst, H_hst);
    // cudaDeviceSynchronize();
    // err = cudaGetLastError();
    // cudaERROR(err);
    // fi = clock();


    st = clock();
    momt_kernelM <<<  L->DimGridMomt_MN, DimBlockMomt_MN, 0, EXECstream[0] >>> (*L);
    momt_kernelN <<<  L->DimGridMomt_MN, DimBlockMomt_MN, 0, EXECstream[1] >>> (*L);
    cudaDeviceSynchronize();
    err = cudaGetLastError();
    cudaERROR(err);
    fi = clock();

    #ifdef DEBUG
        printf("TIME SPENT ON GPU %f\n",(float)(fi-st)/CLOCKS_PER_SEC);
        printf("printing debug information\n" );
        cudaCHK( cudaMemcpy(tmpout, MNout_hst, 2*size_hst[3], cudaMemcpyDeviceToHost) );
        for (size_t i = 0; i < size_hst[2]; i++) {
            if (abs(tmpout[i] - M_f[i]) > ERROR) {
                printf("M[%d,%d] M_cu:%e M_f:%e %e\n", i%size_hst[0], i/size_hst[0], tmpout[i], M_f[i], tmpout[i] - M_f[i]);
            }
        }
        for (size_t i = size_hst[2], j=0; i < 2*size_hst[2]; i++, j++) {
            if (abs(tmpout[i] - N_f[j]) > ERROR) {
                printf("N[%d,%d] N_cu:%e N_f:%e %e\n", (i-size_hst[2])%size_hst[0], (i-size_hst[2])/size_hst[0], tmpout[i], N_f[j], tmpout[i] - N_f[j]);
            }
        }
    #else
        // cudaCHK( cudaMemcpy(M_f, MNout_hst, size_hst[3], cudaMemcpyDeviceToHost) );// FUTURE: delete
        // cudaCHK( cudaMemcpy(N_f, MNout_hst+size_hst[2], size_hst[3], cudaMemcpyDeviceToHost) );// FUTURE: delete
    #endif

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
    const float* __restrict__ R4 = L.R24_hst + L.size_hst[2];
    const float* __restrict__ R5 = L.R35_hst + L.size_hst[1];
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

/*
__global__ void momts_kernel(const float* __restrict__ Z, const float* __restrict__ MN,
                             const float* __restrict__ R24, const float* __restrict__ R35,
                             const float* __restrict__ H) {
    __shared__ float   Z_dat[BLOCKX][BLOCKY+1]; // BLOCKY+1 for preventing bank conflicts
    __shared__ float  MN_dat[BLOCKX][BLOCKY+1][2];
    __shared__ float   H_dat[BLOCKX][BLOCKY+1];

    int row = blockIdx.x*(EXECX) + threadIdx.x;
    int col = blockIdx.y*(EXECY) + threadIdx.y;
    int id = ID(row, col);
    //int idr2 = id+size_dev[2];
    int idmn = id+threadIdx.z*size_dev[2];
    float r1, r2, tot=0.0f, x=0.0f;


    if (row < size_dev[0] && col < size_dev[1]){
    //printf("%d %d %d\n",row, col, id);
        H_dat[threadIdx.x][threadIdx.y] = H[id]; // -->texture mem might be better
        Z_dat[threadIdx.x][threadIdx.y] = Z[id];
        MN_dat[threadIdx.x][threadIdx.y][threadIdx.z] = MN[idmn];
        __syncthreads();


        if (threadIdx.z == 0) { // for M
            if (threadIdx.y != 0 || col == 0) {// boundary condition in threadblocks (if 1st thread is not 1st col)
                if (threadIdx.x < EXECX){ // boundary condition in threadblocks
                    if (row < size_dev[0] - 1){ // IS:IE
                        int ip1 = threadIdx.x+1; // i plus 1
                        int jm1 = threadIdx.y-1; // j minus 1 (only used when col != 0)
                        if (H_dat[threadIdx.x][threadIdx.y] > GX && H_dat[ip1][threadIdx.y] > GX) { // preconditions //NOTE redundant checks
                            r1 = R24[id];
                            r2 = R35[col];
                            tot = MN_dat[threadIdx.x][threadIdx.y][1] + MN_dat[ip1][threadIdx.y][1];
                            x   = -r1*(Z_dat[ip1][threadIdx.y] - Z_dat[threadIdx.x][threadIdx.y]);
                            if (col != 0) {
                                tot += MN_dat[threadIdx.x][jm1][1] + MN_dat[ip1][jm1][1];
                            }else{ // first col
                                tot += tot;
                            }
                            tot *= r2;
                            x += MN_dat[threadIdx.x][threadIdx.y][0];
                            tot += x;
                            if ( tot > EPS || -tot > EPS) MN_out_dev[id] = tot;
                            else MN_out_dev[id] = 0;
                            // if (row == 446 && col == 788) {
                            //     printf("[%d,%d,%d]==================================>momt\n",threadIdx.x,threadIdx.y,threadIdx.z );
                            //     printf("%e\t%e\t%e\t%e\t%e\t%e\n",r1,r2, MN_dat[threadIdx.x][threadIdx.y][1], MN_dat[ip1][threadIdx.y][1], MN_dat[threadIdx.x][jm1][1], MN_dat[ip1][jm1][1] );
                            //     printf("%e\t%e\n",Z_dat[ip1][threadIdx.y], Z_dat[threadIdx.x][threadIdx.y] );
                            // }
                        }
                    }
                }
            }
        }else{ // for N
            if (threadIdx.x != 0 || row == 0) {// boundary condition in threadblocks
                if (threadIdx.y < EXECY){ //boundary condotion in threadblocks
                    if (col < size_dev[1] - 1) { //JS:JE
                        int jp1 = threadIdx.y+1;
                        int im1 = threadIdx.x-1;
                        if (H_dat[threadIdx.x][threadIdx.y] > GX && H_dat[threadIdx.x][jp1] > GX) {  // preconditions
                            r1 = R24[idmn];
                            r2 = R35[col+size_dev[1]];
                            tot  = MN_dat[threadIdx.x][threadIdx.y][0] + MN_dat[threadIdx.x][jp1][0];
                            x    = -r1*(Z_dat[threadIdx.x][jp1] - Z_dat[threadIdx.x][threadIdx.y]);
                            if (row == 0) { // first row
                                tot += tot;
                            }else{
                                tot += MN_dat[im1][threadIdx.y][0] + MN_dat[im1][jp1][0];
                            }
                            tot *= r2;
                            x += MN_dat[threadIdx.x][threadIdx.y][1];
                            tot = x - tot;
                            if ( tot > EPS || -tot > EPS) MN_out_dev[idmn] = tot;
                            else MN_out_dev[idmn] = 0;
                        }
                    }
                }
            }
        }
    }
}
*/