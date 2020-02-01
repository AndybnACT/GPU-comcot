#include "GPUHeader.h"
#include "GPUConfig.h"

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
struct edges_indx{
    int left;
    int right;
};

__global__ void edgeinterp_vert2(const float* __restrict__ outerbase, 
                                 struct edges_indx index, int stride, 
                                 struct GPU_Layer innerL)
{
    const float* __restrict__ edges[2] = {outerbase + index.left, 
                                          outerbase + index.right};
    const uint32_t* __restrict__ size_dev = all_size_dev[innerL.lid];
    const int inner_row[2] = {0, size_dev[0]-1};
    // int row_r = size_dev[0] - 1;
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
        const float * __restrict__ outer = edges[i];
        flux = outer[col*stride];
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



extern "C" void edgeinterp_dbglaunch_(float *mO, float *nO, int *idO, 
                                      float *mA, float *nA, int *idA)
{
    struct GPU_Layer *lO = ldlayer(*idO);
    struct GPU_Layer *lA = ldlayer(*idA);
    struct edges_indx MO_vert;
    float *m_cmp;
    cudaError_t err;

    // edgeinterp_vert <<< lA->DimGrid_JNQ, BLOCKX_JNQ >>> 
    //             (lO->MNdat_hst + lA->corner[0] - 1 - 1, lO->l_size[0], *lA);
    MO_vert.left  = lA->corner[0] - 1 - 1;
    MO_vert.right = lA->corner[1] - 1;
    edgeinterp_vert2 <<< lA->DimGrid_JNQ, BLOCKX_JNQ >>> 
                    (lO->MNdat_hst, MO_vert, lO->l_size[0], *lA);
    cudaDeviceSynchronize();
    err = cudaGetLastError();
    cudaERROR(err);
    
    m_cmp = (float*) malloc(lA->l_size[3]);
    cudaCHK( cudaMemcpy(m_cmp, lA->MNdat_hst, lA->l_size[3], cudaMemcpyDeviceToHost) );
    
    printf("==============lO:%d, lA:%d =============\n", lO->lid, lA->lid);
    for (size_t i = 0; i < lA->l_size[1]; i++) {
        int id = i*lA->l_size[0] + 0;
        if (fabs(m_cmp[id] - mA[id]) >= 1.0e-6) {
            printf("(left) i = %d, %f, %f, diff=%f\n", i, m_cmp[id], mA[id], fabs(m_cmp[id] - mA[id]));
        }
    }
    for (size_t i = 0; i < lA->l_size[1]; i++) {
        int id = i*lA->l_size[0] + lA->l_size[1] - 1;
        // if (m_cmp[id] != 0.0) {
        //     printf("(right !=0, i = %d, %f)\n", i, m_cmp[id]);
        // }
        if (fabs(m_cmp[id] - mA[id]) >= 1.0e-6) {
            printf("(right) i = %d, %f, %f, diff=%f\n", i, m_cmp[id], mA[id], fabs(m_cmp[id] - mA[id]));
        }
    }
    
    free(m_cmp);
}


void jnq(struct GPU_Layer *parent, struct GPU_Layer *child)
{
    
    return;
}


void newq(struct GPU_Layer *parent, struct GPU_Layer *child, int step)
{
    
    return;
}

void jnz(struct GPU_Layer *parent, struct GPU_Layer *child) {
    
    return;
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
