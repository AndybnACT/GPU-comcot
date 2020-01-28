#include "GPUHeader.h"

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
