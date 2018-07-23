#include "GPUHeader.h"
#include "GPUConfig.h"

#include "GPUOpen_BD.h"

extern "C" void openbd_launch_(float *Z_f_complete) {

    cudaError_t err;
    openbd_kernel_tex<<< DimGridOpenBD_LR, DimBlockOpenBD, 0, EXECstream[0] >>>(NULL, H_hst, LEFT);// MN has been changed (:,:,1) <--> (:,:,2)
    openbd_kernel_tex<<< DimGridOpenBD_LR, DimBlockOpenBD, 0, EXECstream[1] >>>(NULL, H_hst, RIGHT);// so use M(:,:,1) directly
    openbd_kernel_tex<<< DimGridOpenBD_TB, DimBlockOpenBD, 0, EXECstream[2] >>>(NULL, H_hst, TOP);
    openbd_kernel_tex<<< DimGridOpenBD_TB, DimBlockOpenBD, 0, EXECstream[3] >>>(NULL, H_hst, BOTTOM);
    cudaDeviceSynchronize();
    err = cudaGetLastError();
    cudaERROR(err);


    #ifdef DEBUG_BOUNDARY_KERNEL
        openbd_kernel<<< DimGridOpenBD_LR, DimBlockOpenBD, 0, EXECstream[0] >>>(MNdat_hst, H_hst, LEFT);// MN has been changed (:,:,1) <--> (:,:,2)
        openbd_kernel<<< DimGridOpenBD_LR, DimBlockOpenBD, 0, EXECstream[1] >>>(MNdat_hst, H_hst, RIGHT);// so use M(:,:,1) directly
        openbd_kernel<<< DimGridOpenBD_TB, DimBlockOpenBD, 0, EXECstream[2] >>>(MNdat_hst, H_hst, TOP);
        openbd_kernel<<< DimGridOpenBD_TB, DimBlockOpenBD, 0, EXECstream[3] >>>(MNdat_hst, H_hst, BOTTOM);
        cudaDeviceSynchronize();
        err = cudaGetLastError();
        cudaERROR(err);
        // printf("printing information for debugging\n" );
        // cudaCHK( cudaMemcpy(tmpout, Zout_hst, size_hst[3], cudaMemcpyDeviceToHost) );
        // for (size_t i = 0; i < size_hst[2]; i++) {
        //     if (abs(tmpout[i] - Z_f_complete[i]) > ERROR) {
        //         printf("Z[%d,%d] Z_cu:%e Z_f:%e %e\n", i%size_hst[0], i/size_hst[0] , tmpout[i], Z_f_complete[i], tmpout[i] - Z_f_complete[i]);
        //     }
        // }
    #endif

}
#ifdef DEBUG_BOUNDARY_KERNEL
    #define ANSCHK(ANS,CHK) if (fabs(ANS-CHK) > ERROR) printf("Z/Z_ inconsistent: Z_%e\tZ:%e\n",ANS, CHK);
#else
    #define ANSCHK(A,B) ;
#endif
__global__ void openbd_kernel(const float* __restrict__ MN, const float* __restrict__ H, bdside BOUNDARY) {
    #define UB 99.0
    float ztmp=0.0, h, m, n, cc;
    uint32_t row=0, col=0;
    switch (BOUNDARY) {
        case RIGHT:
            col = (size_dev[1]-1)*size_dev[0];
            row = blockIdx.x*31*(blockDim.x>>5) + 31*(threadIdx.x>>5) + threadIdx.x%32;
            if (row < size_dev[0]-1) {
                h =  H[row+col];
                m = MN[row+col];// must not load in the following if block, or 1st lane would get 0
                float m_suf = __shfl_up_sync(0xFFFFFFFF,m,1);// must not shuffle in the following if block, or 1st lane would get 0
                n = MN[row+col-size_dev[0]+size_dev[2]];
                cc = 1/sqrtf(GRAV*h);
                if (threadIdx.x % 32 != 0) {
                    if (h > GX) {
                        float uh_2 = 0.25*(m+m_suf)*(m+m_suf);
                        ztmp = sqrtf(n*n + uh_2)*cc;
                        if (n < 0.0) ztmp *= -1;
                        if (ztmp > UB || ztmp < -UB) ztmp = 0.0;
                    }// else {ztmp=0.0;}
                    ANSCHK(Z_out_dev[row+col], ztmp);
                    Z_out_dev[row+col] = ztmp;
                }
                else if (row == 0) { // --|
                    if (h > GX) {
                        ztmp = sqrtf(m*m + n*n)*cc;
                        if (m > 0 || n < 0) ztmp *= -1;
                        if (ztmp > UB || ztmp < -UB) ztmp = 0.0;
                    }
                    ANSCHK(Z_out_dev[col], ztmp);
                    Z_out_dev[col] = ztmp;
                }
            }
            break;
        case LEFT:
            row = blockIdx.x*31*(blockDim.x>>5) + 31*(threadIdx.x>>5) + threadIdx.x%32;
            if (row < size_dev[0]-1) {
                h =  H[row];
                m = MN[row];
                float m_suf = __shfl_up_sync(0xFFFFFFFF,m,1);
                n = MN[row+size_dev[2]];
                cc = 1/sqrtf(GRAV*h);
                if (threadIdx.x % 32 != 0) {
                    if (h > GX) {
                        float uh_2 = 0.25*(m+m_suf)*(m+m_suf);
                        ztmp = sqrtf(n*n + uh_2)*cc;
                        if (n > 0.0) ztmp *= -1;
                        if (ztmp > UB || ztmp < -UB) ztmp = 0.0;
                    }// else {ztmp=0.0;}
                    ANSCHK(Z_out_dev[row], ztmp);
                    Z_out_dev[row] = ztmp;
                }
                else if (row == 0) {//  |--
                    if (h > GX) {
                        ztmp = sqrtf(m*m + n*n)*cc;
                        if (m > 0 || n < 0) ztmp *= -1;
                        if (ztmp > UB || ztmp < -UB) ztmp = 0.0;
                    }
                    ANSCHK(Z_out_dev[col], ztmp);
                    Z_out_dev[col] = ztmp;
                }
            }
            break;
        case TOP: // should use texture
            col = blockIdx.x*31*(blockDim.x>>5) + 31*(threadIdx.x>>5) + threadIdx.x%32;
            if (col < size_dev[1]-1 && col > 0) {
                col *= size_dev[0];
                h =  H[col];
                n = MN[col+size_dev[2]];
                float n_suf = __shfl_up_sync(0xFFFFFFFF,n,1);
                if (threadIdx.x%32 != 0) {
                    if (h > GX) {
                        m = MN[col];
                        cc = 1/sqrtf(GRAV*h);
                        float uh_2 = 0.25*(n+n_suf)*(n+n_suf);

                        ztmp = sqrtf(uh_2 + m*m)*cc;
                        if (m > 0.0) ztmp *= -1;
                        if (ztmp > UB || ztmp < -UB) ztmp = 0.0;
                    }
                    ANSCHK(Z_out_dev[col], ztmp);
                    Z_out_dev[col] = ztmp;
                }
            }
            break;
        case BOTTOM:
            row = size_dev[0]-1;
            col = blockIdx.x*31*(blockDim.x>>5) + 31*(threadIdx.x>>5) + threadIdx.x%32;

            if (col < size_dev[1]-1 && col > 0) { //bottom body
                col *= size_dev[0];
                h =  H[col+row];
                n = MN[col+row+size_dev[2]];
                float n_suf = __shfl_up_sync(0xFFFFFFFF,n,1);
                if (threadIdx.x%32 != 0) {
                    if (h > GX) {
                        m = MN[col+row-1];
                        cc = 1/sqrtf(GRAV*h);
                        float uh_2 = 0.25*(n+n_suf)*(n+n_suf);

                        ztmp = sqrtf(uh_2 + m*m)*cc;
                        if (m < 0.0) ztmp *= -1;
                        if (ztmp > UB || ztmp < -UB) ztmp = 0.0;
                    }
                    ANSCHK(Z_out_dev[col+row], ztmp);
                    Z_out_dev[col+row] = ztmp;
                }
            }else if (col == 0 || col == size_dev[1]-1){ // bottom boundary
                uint32_t id = row+col*size_dev[0];
                h =  H[id];
                if (h > GX) {
                    cc = 1/sqrtf(GRAV*h);
                    m = MN[id-1];
                    if (col == 0) {  // |__
                        n = MN[row+size_dev[2]];
                        ztmp = sqrtf(m*m + n*n);
                        ztmp *= cc;
                        if (m < 0.0 || n > 0.0) ztmp *= -1;
                        if (ztmp > UB || ztmp < -UB) ztmp = 0.0;
                    }else{           //  __|
                        n = MN[row + (col-1)*size_dev[0] + size_dev[2]];
                        ztmp = sqrtf(m*m + n*n);
                        ztmp *= cc;
                        if (m < 0.0 || n < 0.0) ztmp *= -1;
                        if (ztmp > UB || ztmp < -UB) ztmp = 0.0;
                    }
                }
                ANSCHK(Z_out_dev[id], ztmp);
                Z_out_dev[id] = ztmp;
            }
    }
}

__global__ void openbd_kernel_tex(const float* __restrict__ MN, const float* __restrict__ H, bdside BOUNDARY) {
    #define UB 99.0
    float ztmp=0.0, h, m, n, cc;
    uint32_t row=0, col=0;
    switch (BOUNDARY) {
        case RIGHT:
            col = (size_dev[1]-1)*size_dev[0];
            row = blockIdx.x*31*(blockDim.x>>5) + 31*(threadIdx.x>>5) + threadIdx.x%32;
            if (row < size_dev[0]-1) {
                h =  H[row+col];
                m = tex2D<float2>(MNtext, row, size_dev[1]-1).x;//MN[row+col];// must not load in the following if block, or 1st lane would get 0
                float m_suf = __shfl_up_sync(0xFFFFFFFF,m,1);// must not shuffle in the following if block, or 1st lane would get 0
                n = tex2D<float2>(MNtext, row, size_dev[1]-2).y;//MN[row+col-size_dev[0]+size_dev[2]];
                cc = 1/sqrtf(GRAV*h);
                if (threadIdx.x % 32 != 0) {
                    if (h > GX) {
                        float uh_2 = 0.25*(m+m_suf)*(m+m_suf);
                        ztmp = sqrtf(n*n + uh_2)*cc;
                        if (n < 0.0) ztmp *= -1;
                        if (ztmp > UB || ztmp < -UB) ztmp = 0.0;
                    }// else {ztmp=0.0;}
                    Z_out_dev[row+col] = ztmp;
                }
                else if (row == 0) { // --|
                    if (h > GX) {
                        ztmp = sqrtf(m*m + n*n)*cc;
                        if (m > 0 || n < 0) ztmp *= -1;
                        if (ztmp > UB || ztmp < -UB) ztmp = 0.0;
                    }
                    Z_out_dev[col] = ztmp;
                }
            }
            break;
        case LEFT:
            row = blockIdx.x*31*(blockDim.x>>5) + 31*(threadIdx.x>>5) + threadIdx.x%32;
            if (row < size_dev[0]-1) {
                h =  H[row];
                m = tex2D<float2>(MNtext, row, 0).x;//MN[row];
                float m_suf = __shfl_up_sync(0xFFFFFFFF,m,1);
                n = tex2D<float2>(MNtext, row, 0).y;//MN[row+size_dev[2]];
                cc = 1/sqrtf(GRAV*h);
                if (threadIdx.x % 32 != 0) {
                    if (h > GX) {
                        float uh_2 = 0.25*(m+m_suf)*(m+m_suf);
                        ztmp = sqrtf(n*n + uh_2)*cc;
                        if (n > 0.0) ztmp *= -1;
                        if (ztmp > UB || ztmp < -UB) ztmp = 0.0;
                    }// else {ztmp=0.0;}
                    Z_out_dev[row] = ztmp;
                }
                else if (row == 0) {//  |--
                    if (h > GX) {
                        ztmp = sqrtf(m*m + n*n)*cc;
                        if (m > 0 || n < 0) ztmp *= -1;
                        if (ztmp > UB || ztmp < -UB) ztmp = 0.0;
                    }
                    Z_out_dev[col] = ztmp;
                }
            }
            break;
        case TOP: // should use texture
            col = blockIdx.x*31*(blockDim.x>>5) + 31*(threadIdx.x>>5) + threadIdx.x%32;
            if (col < size_dev[1]-1 && col > 0) {
                col *= size_dev[0];
                h =  H[col];
                n = tex2D<float2>(MNtext, 0, col/size_dev[0]).y;//MN[col+size_dev[2]];
                float n_suf = __shfl_up_sync(0xFFFFFFFF,n,1);
                if (threadIdx.x%32 != 0) {
                    if (h > GX) {
                        m = tex2D<float2>(MNtext, 0, col/size_dev[0]).x;//MN[col];
                        cc = 1/sqrtf(GRAV*h);
                        float uh_2 = 0.25*(n+n_suf)*(n+n_suf);

                        ztmp = sqrtf(uh_2 + m*m)*cc;
                        if (m > 0.0) ztmp *= -1;
                        if (ztmp > UB || ztmp < -UB) ztmp = 0.0;
                    }
                    Z_out_dev[col] = ztmp;
                }
            }
            break;
        case BOTTOM:
            row = size_dev[0]-1;
            col = blockIdx.x*31*(blockDim.x>>5) + 31*(threadIdx.x>>5) + threadIdx.x%32;

            if (col < size_dev[1]-1 && col > 0) { //bottom body
                col *= size_dev[0];
                h =  H[col+row];
                n = tex2D<float2>(MNtext, row, col/size_dev[0]).y;//MN[col+row+size_dev[2]];
                float n_suf = __shfl_up_sync(0xFFFFFFFF,n,1);
                if (threadIdx.x%32 != 0) {
                    if (h > GX) {
                        m = tex2D<float2>(MNtext, row-1, col/size_dev[0]).x;//MN[col+row-1];
                        cc = 1/sqrtf(GRAV*h);
                        float uh_2 = 0.25*(n+n_suf)*(n+n_suf);

                        ztmp = sqrtf(uh_2 + m*m)*cc;
                        if (m < 0.0) ztmp *= -1;
                        if (ztmp > UB || ztmp < -UB) ztmp = 0.0;
                    }
                    Z_out_dev[col+row] = ztmp;
                }
            }else if (col == 0 || col == size_dev[1]-1){ // bottom boundary
                uint32_t id = row+col*size_dev[0];
                h =  H[id];
                if (h > GX) {
                    cc = 1/sqrtf(GRAV*h);
                    m = tex2D<float2>(MNtext, row-1, col).x;//MN[id-1];
                    if (col == 0) {  // |__
                        n = tex2D<float2>(MNtext, row, 0).y;//MN[row+size_dev[2]];
                        ztmp = sqrtf(m*m + n*n);
                        ztmp *= cc;
                        if (m < 0.0 || n > 0.0) ztmp *= -1;
                        if (ztmp > UB || ztmp < -UB) ztmp = 0.0;
                    }else{           //  __|
                        n = tex2D<float2>(MNtext, row, col-1).y;//MN[row + (col-1)*size_dev[0] + size_dev[2]];
                        ztmp = sqrtf(m*m + n*n);
                        ztmp *= cc;
                        if (m < 0.0 || n < 0.0) ztmp *= -1;
                        if (ztmp > UB || ztmp < -UB) ztmp = 0.0;
                    }
                }
                Z_out_dev[id] = ztmp;
            }
    }
}
