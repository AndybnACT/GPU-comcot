#include "GPUHeader.h"
#include "GPUConfig.h"
#include "debug_option.h"

void show_grid(struct GPU_Layer *L)
{
    printf("lid = %d\n", L->lid );
    printf("parent = %d\n", L->plid );
    printf("level = %d\n", L->lvl);
    printf("size = %d x %d\n", L->l_size[0], L->l_size[1]);
}

extern "C" void for_alll_grids_rcsv(void(*__do)(struct GPU_Layer *), struct GPU_Layer *l)
{
    if (!l) {
        return;
    }
    
    __do(l);
    for_alll_grids_rcsv(__do, l->child);
    for_alll_grids_rcsv(__do, l->sibling);
    
    return;
}

extern "C" void gcomcot_show_grids_()
{
    for_alll_grids_rcsv(show_grid, Root_Grid);
    return;
}

// nx row
// ny col
__global__ void edgeinterp_vert(const float* __restrict__ outer, int stride, 
                                struct GPU_Layer innerL)
{
    const uint32_t __restrict__ *size_dev = all_size_dev[innerL.lid];
    int inner_row = 0;
    // int row_r = size_dev[0] - 1;
    int col = threadIdx.x%32 + 30*(threadIdx.x>>5) + blockIdx.x*30*(blockDim.x>>5) - 1;
    int col_end = innerL.corner[3];
    const int ir = innerL.rel_size;
    
    int inner_col = col*ir + 1;
    
    col += innerL.corner[2] - 1;
    
    float flux;
    double flux_s, flux_e, c1, c2, nn;
    
    nn = 2.0 * (double) ir;
    
    flux = outer[col*stride + innerL.corner[0] - 1 - 1];
    flux_e = (double) __shfl_down_sync(0xFFFFFFFF, flux, 1);
    flux_s = (double) __shfl_up_sync(0xFFFFFFFF, flux, 1);
    flux_e = 0.5*(flux_e + (double) flux);
    flux_s = 0.5*(flux_s + (double) flux);
    
    if (threadIdx.x%32 == 0 || threadIdx.x%32 == 31) {
        return;
    }
    
    if (col < col_end) {
        c1 = (flux_e-flux_s)/nn;
        c2 = -c1 + flux_s;
        
        for (size_t i = 0; i < ir; i++) {
            innerL.MNdat_hst[ID(inner_row, inner_col)] = 2.0*((float)(i+1))*c1 + c2;
            if (innerL.H_hst[ID(inner_row, inner_col)] + innerL.Zdat_hst[ID(inner_row, inner_col)] <= GX)
                innerL.MNdat_hst[ID(inner_row, inner_col)] = 0.0;
            inner_col++;
        }
        
    }
}

__device__ float newq_vert(struct GPU_Layer outerL, struct GPU_Layer innerL,
                          int row, int col, int t, int side)
{
    const uint32_t* __restrict__ size_dev = all_size_dev[outerL.lid];
    float c1, c2;
    float n, n_ip1, n_jm1, n_ip1_jm1, tot_n;
    float m, xm;
    float z, z_ip1, h, h_ip1;
    float r2, r3;
    float yflux = 0.0f;

    c1 = ((float) t - 1.0f)/((float) innerL.rel_time);
    c2 = 1.0 - c1;

    h = outerL.H_hst[ID(row, col)];
    h_ip1 =  outerL.H_hst[ID(row + 1, col)];
    z = outerL.Zout_hst[ID(row, col)];
    z_ip1 = outerL.Zout_hst[ID(row + 1, col)];

    if (t == 2) {
        n = outerL.MNdat_hst[ID2E(row, col, 1)];
        n_ip1 = outerL.MNdat_hst[ID2E(row + 1, col, 1)];
        n_jm1 = __shfl_up_sync(0xFFFFFFFF, n, 1);
        n_ip1_jm1 = __shfl_up_sync(0xFFFFFFFF, n_ip1, 1);
        if (threadIdx.x%32 == 0) { // first lane cannot get jm1 via shfl_up
            n_jm1 = outerL.MNdat_hst[ID2E(row, col - 1, 1)];
            n_ip1_jm1 = outerL.MNdat_hst[ID2E(row + 1, col - 1, 1)];
        }
    }

    if (col < innerL.corner[2] - 2 || col > innerL.corner[3] )
        return 0;

    if (h + z > GX && h_ip1 + z_ip1 > GX) {
        m = outerL.MNdat_hst[ID2E(row, col, 0)];
        if (t == 2) {
            r2 = outerL.R24_hst[ID(row, col)];
            r3 = outerL.R35_hst[col];

            tot_n = n + n_ip1 + n_jm1 + n_ip1_jm1;
            xm = m - r2*(z_ip1 - z) + r3*tot_n;
            if (fabs(xm) < EPS)
                xm = 0.0f;

            outerL.yflux[col + outerL.l_size[1]*side] = xm;
            if (col == 1411 && side == 0) {
                printf("%d, %d, %e, %e, %d, %d\n", outerL.lid, innerL.lid, xm, tot_n, threadIdx.x, blockIdx.x);
            }
        }
        yflux = c1*outerL.yflux[col + outerL.l_size[1]*side] + c2*m;
    }

    return yflux;
}

__device__ float newq_hori(struct GPU_Layer outerL, struct GPU_Layer innerL,
                          int row, int col, int t, int side)
{
    const uint32_t* __restrict__ size_dev = all_size_dev[outerL.lid];
    float c1, c2;
    float m, m_im1, m_jp1, m_im1_jp1, tot_m;
    float n, xn;
    float z, z_jp1, h, h_jp1;
    float r4, r5;
    float xflux = 0.0f;

    c1 = ((float) t - 1.0f)/((float) innerL.rel_time);
    c2 = 1.0 - c1;

    h = outerL.H_hst[ID(row, col)];
    h_jp1 =  outerL.H_hst[ID(row, col + 1)];
    z = outerL.Zout_hst[ID(row, col)];
    z_jp1 = outerL.Zout_hst[ID(row, col + 1)];

    if (t == 2) {
        m = outerL.MNdat_hst[ID(row, col)];
        m_jp1 = outerL.MNdat_hst[ID(row, col + 1)];
        m_im1 = __shfl_up_sync(0xFFFFFFFF, m, 1);
        m_im1_jp1 = __shfl_up_sync(0xFFFFFFFF, m_jp1, 1);
        if (threadIdx.x%32 == 0) {
            m_im1 = outerL.MNdat_hst[ID(row - 1, col)];
            m_im1_jp1 = outerL.MNdat_hst[ID(row - 1, col + 1)];
        }
    }

    if (row < innerL.corner[0] - 2 || row > innerL.corner[1])
        return 0;

    if (h + z > GX && h_jp1 + z_jp1 > GX) {
        n = outerL.MNdat_hst[ID2E(row, col, 1)];
        if (t == 2) {
            r4 = outerL.R24_hst[ID2E(row, col, 1)];
            r5 = outerL.R35_hst[col + size_dev[1]];

            tot_m = m_im1 + m_im1_jp1 + m + m_jp1;
            xn = n - r4*(z_jp1 - z) - r5*tot_m;
            if (fabs(xn) < EPS)
                xn = 0.0f;

            outerL.xflux[row + outerL.l_size[0]*side] = xn;
        }
        xflux = c1*outerL.xflux[row + outerL.l_size[0]*side] + c2*n;
    }

    return xflux;
}

__global__ void edgeinterp_vert2(struct GPU_Layer innerL, struct GPU_Layer outerL, int T)
{
    const uint32_t* __restrict__ size_dev = all_size_dev[innerL.lid];
    const int outer_row[2] = {innerL.corner[0] - 1 - 1,
                              innerL.corner[1] - 1};
    const int inner_row[2] = {0, size_dev[0]-1};
    const int colid = threadIdx.x%32 + 30*(threadIdx.x>>5) + blockIdx.x*30*(blockDim.x>>5) - 1;
    int col = colid;
    int col_end = innerL.corner[3];
    const int ir = innerL.rel_size;
    
    col += innerL.corner[2] - 1;
    
    float flux;
    float flux_s, flux_e, c1[2], c2[2], nn;
    
    nn = 2.0 * (float) ir;
    
#pragma unroll 2
    for (size_t i = 0; i < 2; i++) {
        if (T == 1) {
            flux = outerL.MNdat_hst[outer_row[i] + col*outerL.l_size[0]];
        }
        else {
            flux = newq_vert(outerL, innerL, outer_row[i], col, T, i);
            __syncwarp();
        }
        flux_e = __shfl_down_sync(0xFFFFFFFF, flux, 1);
        flux_s = __shfl_up_sync(0xFFFFFFFF, flux, 1);
        flux_e = 0.5*(flux_e + flux);
        flux_s = 0.5*(flux_s + flux);
        c1[i] = (flux_e-flux_s)/nn;
        c2[i] = -c1[i] + flux_s;
        // if (colid == 59 || colid == 89 || colid == 209 || colid == 269 || colid == 29) {
        //     printf("row = %d col = %d colid = %d, threadIdx.x=%d, fluxes=%f, %f, %f\n",
        //             outer_row[i], col, colid, threadIdx.x, flux, flux_s, flux_e);
        // }
    }

    // The first and last lanes cannot get flux via shfl_up/down_sync
    if (threadIdx.x%32 == 0 || threadIdx.x%32 == 31) {
        return;
    }
    
    if (col < col_end) {
#pragma unroll 2
        for (size_t i = 0; i < 2; i++) {
            int inner_col = colid*ir + 1;
            for (size_t j = 0; j < ir; j++) {
                innerL.MNdat_hst[ID(inner_row[i], inner_col)] = 
                    2.0*((float)(j+1))*c1[i] + c2[i];
                if (innerL.H_hst[ID(inner_row[i], inner_col)] + 
                    innerL.Zdat_hst[ID(inner_row[i], inner_col)] <= GX)
                        innerL.MNdat_hst[ID(inner_row[i], inner_col)] = 0.0;
                inner_col++;
            }
        }
    }
}

__global__ void edgeinterp_hori2(struct GPU_Layer innerL, struct GPU_Layer outerL, int T)
{
    const uint32_t* __restrict__ size_dev = all_size_dev[innerL.lid];
    const int outer_col[2] = {(innerL.corner[2] - 1 - 1),
                              (innerL.corner[3] - 1)};
    const int inner_col[2] = {0, size_dev[1]-1};
    const int rowid = threadIdx.x%32 + 30*(threadIdx.x>>5) + blockIdx.x*30*(blockDim.x>>5) - 1;
    int row = rowid;
    int row_end = innerL.corner[1];
    const int ir = innerL.rel_size;
    
    row += innerL.corner[0] - 1;
    
    float flux;
    float flux_s, flux_e, c1[2], c2[2], nn;
    
    nn = 2.0 * (float) ir;
    
#pragma unroll 2
    for (size_t i = 0; i < 2; i++) {
        if (T == 1) {
            flux = outerL.MNdat_hst[outerL.l_size[2] + outer_col[i]*outerL.l_size[0] + row];
        }
        else {
            flux = newq_hori(outerL, innerL, row, outer_col[i], T, i);
            __syncwarp();
        }
        flux_e = __shfl_down_sync(0xFFFFFFFF, flux, 1);
        flux_s = __shfl_up_sync(0xFFFFFFFF, flux, 1);
        flux_e = 0.5*(flux_e + flux);
        flux_s = 0.5*(flux_s + flux);
        c1[i] = (flux_e-flux_s)/nn;
        c2[i] = -c1[i] + flux_s;
    }
    
    if (threadIdx.x%32 == 0 || threadIdx.x%32 == 31) {
        return;
    }
    
    if (row < row_end) {
#pragma unroll 2
        for (size_t i = 0; i < 2; i++) {
            int inner_row = rowid*ir + 1;
            for (size_t j = 0; j < ir; j++) {
                innerL.MNdat_hst[ID2E(inner_row, inner_col[i], 1)] = 
                    2.0*((float)(j+1))*c1[i] + c2[i];
                if (innerL.H_hst[ID(inner_row, inner_col[i])] + 
                    innerL.Zdat_hst[ID(inner_row, inner_col[i])] <= GX)
                        innerL.MNdat_hst[ID2E(inner_row, inner_col[i], 1)] = 0.0;
                inner_row++;
            }
        }
    }    
}

extern "C" void edgeinterp_dbglaunch_(float *mO, float *nO, int *idO, 
                                      float *mA, float *nA, int *idA, int *t,
                                      float *yflux, float *xflux)
{
    clock_t st, fi;
    struct GPU_Layer *lO = ldlayer(*idO);
    struct GPU_Layer *lA = ldlayer(*idA);
    float *mn_cmp, *yflux_cmp, *xflux_cmp;
    cudaError_t err;

    st = clock();
    // edgeinterp_vert <<< lA->DimGrid_JNQ, BLOCKX_JNQ >>> 
    //             (lO->MNdat_hst + lA->corner[0] - 1 - 1, lO->l_size[0], *lA);
    edgeinterp_vert2 <<< lA->DimGrid_JNQV, BLOCKX_JNQ, 0, EXECstream[0] >>> 
                    (*lA, *lO, *t);
    
    edgeinterp_hori2 <<< lA->DimGrid_JNQH, BLOCKX_JNQ, 0, EXECstream[1] >>>
                    (*lA, *lO, *t);
    
    cudaDeviceSynchronize();
    err = cudaGetLastError();
    cudaERROR(err);
    fi = clock();
    printf("TIME SPENT ON GPU (EDGE_INTERP->%d) %f\n",(float)(fi-st)/CLOCKS_PER_SEC, lA->lid);
    
    mn_cmp = (float*) malloc(2*lA->l_size[3]);
    xflux_cmp = (float*) malloc(2*lO->l_size[0]*sizeof(float));
    yflux_cmp = (float*) malloc(2*lO->l_size[1]*sizeof(float));
    cudaCHK( cudaMemcpy(mn_cmp, lA->MNdat_hst, 2*lA->l_size[3], cudaMemcpyDeviceToHost) );
    cudaCHK( cudaMemcpy(yflux_cmp, lO->yflux, 2*lO->l_size[1]*sizeof(float), cudaMemcpyDeviceToHost) );
    cudaCHK( cudaMemcpy(xflux_cmp, lO->xflux, 2*lO->l_size[0]*sizeof(float), cudaMemcpyDeviceToHost) );
    
    printf("==============lO:%d, lA:%d =============\n", lO->lid, lA->lid);
    if (*t == 2) {
        for (size_t i = lA->corner[2]-1; i <= lA->corner[3]; i++) {
            int id = i;
            if (assert_diff(yflux_cmp[id], yflux[id])) {
                printf("(left-flux) i = %d, %f, %f, diff=%f\n", i, yflux_cmp[id], yflux[id], fabs(yflux_cmp[id] - yflux[id]));
            }
        }

        for (size_t i = lA->corner[2]-1; i <= lA->corner[3]; i++) {
            int id = i + lO->l_size[1];
            if (assert_diff(yflux_cmp[id], yflux[id])) {
                printf("(right-flux) i = %d, %f, %f, diff=%f\n", i, yflux_cmp[id], yflux[id], fabs(yflux_cmp[id] - yflux[id]));
            }
        }

        for (size_t i = lA->corner[0]-1; i < lA->corner[1]; i++) {
            int id = i;
            if (assert_diff(xflux_cmp[id], xflux[id])) {
                printf("(bottom-flux) i = %d, %f, %f, diff=%f\n", i, xflux_cmp[id], xflux[id], fabs(xflux_cmp[id] - xflux[id]));
            }
        }

        for (size_t i = lA->corner[0]-1; i <= lA->corner[1]; i++) {
            int id = i + lO->l_size[0];
            if (assert_diff(xflux_cmp[id], xflux[id])) {
                printf("(top-flux) i = %d, %f, %f, diff=%f\n", i, xflux_cmp[id], xflux[id], fabs(xflux_cmp[id] - xflux[id]));
            }
        }
    }

    for (size_t i = 0; i < lA->l_size[1]; i++) {
        int id = i*lA->l_size[0] + 0;
        if (assert_diff(mn_cmp[id], mA[id])) {
            printf("(left) i = %d, %f, %f, diff=%f\n", i, mn_cmp[id], mA[id], fabs(mn_cmp[id] - mA[id]));
        }
    }
    for (size_t i = 0; i < lA->l_size[1]; i++) {
        int id = i*lA->l_size[0] + lA->l_size[0] - 1;
        if (assert_diff(mn_cmp[id], mA[id]) ) {
            printf("(right) i = %d, %f, %f, diff=%f\n", i, mn_cmp[id], mA[id], fabs(mn_cmp[id] - mA[id]));
        }
    }
    
    for (size_t i = 0; i < lA->l_size[0]; i++) {
        int id = i;
        if (assert_diff(mn_cmp[id + lA->l_size[2]], nA[id])) {
            printf("(bottom) i = %d, %f, %f, diff=%f\n", i, mn_cmp[id + lA->l_size[2]], 
                    nA[id], fabs(mn_cmp[id + lA->l_size[2]] - nA[id]));
        }
    }
    
    for (size_t i = 0; i < lA->l_size[0]; i++) {
        int id = i + lA->l_size[0]*(lA->l_size[1] - 1);
        if (assert_diff(mn_cmp[id + lA->l_size[2]], nA[id])) {
            printf("(top) i = %d, %f, %f, diff=%f\n", i, mn_cmp[id + lA->l_size[2]], 
                    nA[id], fabs(mn_cmp[id + lA->l_size[2]] - nA[id]));
        }
    }
    
    free(mn_cmp);
    free(yflux_cmp);
    free(xflux_cmp);
}

void jnq(struct GPU_Layer *lO, struct GPU_Layer *lA)
{
    clock_t st, fi;
    cudaError_t err;

    st = clock();
    // edgeinterp_vert <<< lA->DimGrid_JNQ, BLOCKX_JNQ >>> 
    //             (lO->MNdat_hst + lA->corner[0] - 1 - 1, lO->l_size[0], *lA);
    edgeinterp_vert2 <<< lA->DimGrid_JNQV, BLOCKX_JNQ, 0, EXECstream[0] >>> 
                    (*lA, *lO, 1);
    
    edgeinterp_hori2 <<< lA->DimGrid_JNQH, BLOCKX_JNQ, 0, EXECstream[1] >>>
                    (*lA, *lO, 1);
    
    cudaDeviceSynchronize();
    err = cudaGetLastError();
    cudaERROR(err);
    fi = clock();
    return;
}


void newq(struct GPU_Layer *lO, struct GPU_Layer *lA, int step)
{
    clock_t st, fi;
    cudaError_t err;

    st = clock();
    // edgeinterp_vert <<< lA->DimGrid_JNQ, BLOCKX_JNQ >>> 
    //             (lO->MNdat_hst + lA->corner[0] - 1 - 1, lO->l_size[0], *lA);
    edgeinterp_vert2 <<< lA->DimGrid_JNQV, BLOCKX_JNQ, 0, EXECstream[0] >>> 
                    (*lA, *lO, step);
    
    edgeinterp_hori2 <<< lA->DimGrid_JNQH, BLOCKX_JNQ, 0, EXECstream[1] >>>
                    (*lA, *lO, step);
    
    cudaDeviceSynchronize();
    err = cudaGetLastError();
    cudaERROR(err);
    fi = clock();
    return;
}

__global__ void jnz_kernel_deprecated(float * __restrict__ Z_out, int stride,
                           float * __restrict__ Z_in, int nr, int nc,
                           int rel_size, int thread_xbound)
{
    __shared__ float col_reduct[32];
    int row, col, col_end;
    float row_reduct, final_val;
    
    // row = threadIdx.x + blockIdx.x*(32/rel_size*rel_size) ;
    row = threadIdx.x + blockIdx.x*thread_xbound;
    col = threadIdx.y + blockIdx.y*10*rel_size;
    col_end = col + 10 < nc ? col + 10 : nc;
    
    for (size_t i = 0, j = col; j < col_end; i++, j++) {
        if (row < nr) {
            row_reduct = Z_in[j*nr + row];
        }
        if (rel_size == 2) {
            row_reduct += __shfl_down_sync(0xFFFFFFFF, row_reduct, 1);
            if (threadIdx.x & 0x1 == 0) { // x%2 == 0
                col_reduct[threadIdx.y*16 + threadIdx.x/2] = row_reduct;
            }
            __syncthreads();
            if (threadIdx.x & 0x1 == 0 && threadIdx.y == 0) {
                final_val = row_reduct;
                final_val += col_reduct[threadIdx.x/2 + 16];
                Z_out[stride*i + threadIdx.x/2] = final_val/4.0f;
            }
        }
        if (rel_size == 4) {
            row_reduct += __shfl_down_sync(0xFFFFFFFF, row_reduct, 2);
            row_reduct += __shfl_down_sync(0xFFFFFFFF, row_reduct, 1);
            if (threadIdx.x & 0x3 == 0) { // x%4 == 0
                /* code */
            }
        }
        if (rel_size == 8) {
            row_reduct += __shfl_down_sync(0xFFFFFFFF, row_reduct, 4);
            row_reduct += __shfl_down_sync(0xFFFFFFFF, row_reduct, 2);
            row_reduct += __shfl_down_sync(0xFFFFFFFF, row_reduct, 1);
            if (threadIdx.x & 0x7 == 0) { // x%8 == 0
                
            }
        }
    }
}

#define WRP_ALL_MASK 0xffffffff
template <int relsize>
__global__ void jnz_wrap_based_kernel(const struct GPU_Layer LO,
                                      const struct GPU_Layer LA,
                                      const int reltime)
{
    const float half = relsize * relsize / 2.0;
    float* __restrict__ ZO = LO.Zout_hst;
    const float* __restrict__ ZA1 = LA.Zdat_hst;
    const float* __restrict__ ZA2 = LA.Zout_hst;
    const float* __restrict__ HA = LA.H_hst;
    int row, col, row_end, col_end;
    int dest_row, dest_col;
    const int nr_row_lo = LO.l_size[0];
    const int nr_row_la = LA.l_size[0];
    const int offset_i = LA.corner[0] - 1; // translate to C's zero-based index
    const int offset_j = LA.corner[2] - 1;
    bool writer;
    float z, h;
    int cnt;
    float res;
    float reduct_z;
    int reduct_cnt;

    row = threadIdx.x + blockIdx.x * blockDim.x + 1;
    col = blockIdx.y * LOOPDEPTH + 1; // threadIdx.y is always 0
    row_end = LA.l_size[0];
    col_end = (blockIdx.y == gridDim.y - 1) ? LA.l_size[1]: 
                                              col + LOOPDEPTH;
    if (row >= row_end)
        return;
    
    dest_row = offset_i + (threadIdx.x + blockIdx.x * blockDim.x) / relsize;
    dest_col = offset_j + (blockIdx.y * LOOPDEPTH) / relsize;
    writer = threadIdx.x % relsize == 0;

    while (col < col_end) {
        int c = col;
        reduct_z = 0;
        reduct_cnt = 0;
        for (; c < col + relsize; c++) {
            z = ZA2[ID_MATRIX(row, c, nr_row_la)];
            h = HA[ID_MATRIX(row, c, nr_row_la)]; 
            if (z + h <= GX) {
                z = 0;
                cnt = 0;
            } else {
                if (reltime % 2 == 0) 
                    z = (z + ZA1[ID_MATRIX(row, c, nr_row_la)]) / 2.0; 
                cnt = 1;
            }
            /*
            if (dest_row == 391 && dest_col == 1299) {
                printf("z = %e, h = %e, row = %d, col = %d, cnt = %d, tidx = %d, bidy = %d\n",
                        z, h, row, c, cnt, threadIdx.x, blockIdx.y);
            }
            */
            // reduction starts
            __syncwarp();
            if (relsize == 32) {
                z += __shfl_down_sync(WRP_ALL_MASK, z, 16);
                cnt += __shfl_down_sync(WRP_ALL_MASK, cnt, 16);
            }
            if (relsize >= 16) {
                z += __shfl_down_sync(WRP_ALL_MASK, z, 8);
                cnt += __shfl_down_sync(WRP_ALL_MASK, cnt, 8);
            }
            if (relsize >= 8) {
                z += __shfl_down_sync(WRP_ALL_MASK, z, 4);
                cnt += __shfl_down_sync(WRP_ALL_MASK, cnt, 4);
            }
            if (relsize >= 4) {
                z += __shfl_down_sync(WRP_ALL_MASK, z, 2);
                cnt += __shfl_down_sync(WRP_ALL_MASK, cnt, 2);
            }
            if (relsize >= 2) {
                z += __shfl_down_sync(WRP_ALL_MASK, z, 1);
                cnt += __shfl_down_sync(WRP_ALL_MASK, cnt, 1);
            }
            reduct_z += z;
            reduct_cnt += cnt;
        }
        // write back
        if (writer) {
            /*
            if (dest_row == 391 && dest_col == 1299) {
                printf("reducted_z = %e, reducted_cnt = %d, %d/%d, %d/%d\n",
                        reduct_z, reduct_cnt, 
                        threadIdx.x, blockIdx.x, threadIdx.y, blockIdx.y);
            }
            */
            res = (reduct_cnt > half) ? (reduct_z / (float) reduct_cnt) : 0.0;
            ZO[ID_MATRIX(dest_row, dest_col, nr_row_lo)] = res;
            dest_col++;
        }
        col = c; // c would be
    }

    return;
}

void jnz(struct GPU_Layer *L_p, struct GPU_Layer *L_c) 
{
    clock_t st, fi;
    cudaError_t err;
    
    st = clock();
    switch (L_c->rel_size) {
    case 2:
        jnz_wrap_based_kernel<2> <<< L_c->DimGridJnz, L_c->DimBlockJnz >>>
                                    (*L_p, *L_c, L_c->rel_time);  
        break;
    case 4: 
        jnz_wrap_based_kernel<4> <<< L_c->DimGridJnz, L_c->DimBlockJnz >>>
                                    (*L_p, *L_c, L_c->rel_time);  
        break;
    case 8:  
        jnz_wrap_based_kernel<8> <<< L_c->DimGridJnz, L_c->DimBlockJnz >>>
                                    (*L_p, *L_c, L_c->rel_time);  
        break;
    case 16:  
        jnz_wrap_based_kernel<16> <<< L_c->DimGridJnz, L_c->DimBlockJnz >>>
                                    (*L_p, *L_c, L_c->rel_time);  
        break;
    case 32: 
        jnz_wrap_based_kernel<32> <<< L_c->DimGridJnz, L_c->DimBlockJnz >>>
                                    (*L_p, *L_c, L_c->rel_time);  
        break;
    } 
    cudaDeviceSynchronize();
    err = cudaGetLastError();
    cudaERROR(err);
    fi = clock();
    
#ifdef DEBUG_CORE
    printf("TIME SPENT ON GPU (JNZ) %f\n",(float)(fi-st)/CLOCKS_PER_SEC);
#endif /* DEBUG_CORE */
    return;
}

extern "C" void jnz_dbglaunch_(float *Z_f, int *outid, int *inid) {
    GPU_Layer *parent = ldlayer(*outid);
    GPU_Layer *child  = ldlayer(*inid);
    float *z_cmp = (float*) malloc(parent->l_size[3]);
    float *z_ori = (float*) malloc(child->l_size[3]);
    
    jnz(parent, child);
    cudaCHK( cudaMemcpy(z_cmp, parent->Zout_hst, parent->l_size[3], cudaMemcpyDeviceToHost) );
    cudaCHK( cudaMemcpy(z_ori, child->Zout_hst, child->l_size[3], cudaMemcpyDeviceToHost) );
    printf("child corner %d %d %d %d\n", 
        child->corner[0], child->corner[1], child->corner[2], child->corner[3]);
    CMP_VAR(z_cmp, Z_f, 
            parent->l_size[0], parent->l_size[1],
            "JNZ, layerid=%d", child->lid);
    
    free(z_ori);
    free(z_cmp);
}

void update(struct GPU_Layer *L)
{
    
    // maybe we could place cudaDeviceSynchronize() before some kernel launch 
    // so that we can run C/Fortran code simutaneously with cuda code
    cuda_update_layer_(&L->lid);
    
    return;
}

void subgrid_solver_rcsv(struct GPU_Layer *current)
{
    struct GPU_Layer *parent, *sibling, *child;
    const int *id = &current->lid;
    
    int halfway = 0;
    
    if (!current) {
        return;
    }
    
    if (current->plid == -1) {
        return;
    }
    parent = ldlayer(current->plid);
    
    sibling = current->sibling;
    child = current->child;
    
    halfway = (current->rel_time/2) + 1;
    for (size_t i = 1; i <= current->rel_time; i++) {
        if (i == 1) {
            jnq(parent, current);
        }
        else {
            newq(parent, current, i);
        }
        
        mass_launch_(NULL , NULL, NULL, id);
        
        subgrid_solver_rcsv(current->child);

        momt_launch_(NULL, NULL, NULL, id);
        
        if (i == halfway) {
            jnz(parent, current);
        }
        
        update(current);
    }
    
    subgrid_solver_rcsv(sibling);
    
    return;
}

extern "C" void all_grid_launch_(void)
{
    // All subgrids start from 2. 1 is the Root_Grid.
    return subgrid_solver_rcsv(ldlayer(2));
}

