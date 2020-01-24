#include "GPUHeader.h"
#include "GPUConfig.h"

// global variables-----------------------------
struct GPU_Layer Layer_struct[MAX_LAYERS];
uint32_t all_size[MAX_LAYERS][4];
__constant__ __device__ uint32_t all_size_dev[MAX_LAYERS][4];

float *Zout_hst, *MNout_hst;
float *MNdat_hst, *Zdat_hst;
float *R24_hst, *R35_hst, *H_hst;
float *R_MASS_hst;
float *Zmax_hst;
// __device__ float *R35_dev;
// __device__ float *R24_dev, *H_dev;
// __device__ float *Z_dat_dev, *MN_dat_dev;
// texture<float, cudaTextureType2D, cudaReadModeElementType> ZtexRef;
// texture<float, cudaTextureType2D, cudaReadModeElementType> MNtexRef;


cudaStream_t EXECstream[NUMSTREAM];
dim3 DimBlockMomt_MN(BLOCKX_MOMT,1,1);
dim3 DimGridMomt_MN(0,1,1);
dim3 DimBlockMomt(BLOCKX,BLOCKY,2);
dim3 DimGridMomt(0,0,1);
dim3 DimBlockMass(BLOCKX_MASS,1,1);
dim3 DimGridMass(1,1,1);
dim3 DimBlockOpenBD(BLOCKX_OPENBD,1,1);
dim3 DimGridOpenBD_LR(0,1,1);
dim3 DimGridOpenBD_TB(0,1,1);
size_t GridMaxAmp;

cudaDeviceProp dev_prop;
#ifdef DEBUG
    float *tmpout;
#endif
//----------------------------------------------
extern "C" void cuda_update_layer_(int *);
void cudaMalloc2E(void**, void*, void*, size_t);
#ifdef DEBUG
    __host__ __device__ void prt_mat(float*, size_t, size_t);
    __global__ void GPUmemtest(void);
    __global__ void CHECK_VAR(void);
    extern "C" void cmpvar_(const float *, int *);
#endif

extern "C" void cuda_alloc_layer(struct GPU_Layer *L, float *R1_f, float *R2_f, 
                                 float *R3_f, float *R4_f, float *R5_f, float *R6_f, 
                                 float *R11_f, float *H_f, float *Z_f, int *row, int *col){

    //float *R24_hst, *R35_hst, *H_hst;
    float *R3, *R5;
    float *R_MASS;



    L->size_hst[0] = *row;
    L->size_hst[1] = *col;
    L->size_hst[2] = (*row)*(*col);
    L->size_hst[3] = (*row)*(*col)*sizeof(float);
    // allocate variables
    // const === R2 R4 H ===
    cudaCHK( cudaMalloc(&L->H_hst, L->size_hst[3]) );

    cudaMalloc2E((void**)&L->R24_hst, R2_f, R4_f, L->size_hst[3]);
    //cudaCHK( cudaMemcpyToSymbol(R24_dev, &R24_hst, sizeof(float*)) );
    // const variables R3 R5, R1 R6 R11
    size_t R_size = sizeof(float)*L->size_hst[1];
    //R1 = (float*) malloc(R_size);
    R3 = (float*) malloc(R_size);
    R5 = (float*) malloc(R_size);
    //R11 =(float*) malloc(R_size);
    R_MASS = (float*) malloc(4*R_size);

    for (size_t i = 0; i < L->size_hst[1]; i++) {
        R3[i] = R3_f[i*L->size_hst[0]];
        R5[i] = R5_f[i*L->size_hst[0]];
        R_MASS[4*i]   = R1_f[i*L->size_hst[0]];
        R_MASS[4*i+1] = R6_f[i];
        R_MASS[4*i+2] = R11_f[i*L->size_hst[0]];
    }


    // cudaCHK( cudaMalloc(&R1_hst, R_size) );
    // cudaCHK( cudaMemcpy(R1_hst, R1, R_size, cudaMemcpyHostToDevice) );

    // cudaCHK( cudaMalloc(&R6_hst, R_size) );
    // cudaCHK( cudaMemcpy(R6_hst, R6_f, R_size, cudaMemcpyHostToDevice) );

    // cudaCHK( cudaMalloc(&R11_hst, R_size) );
    // cudaCHK( cudaMemcpy(R11_hst, R11, R_size, cudaMemcpyHostToDevice) );
    cudaCHK( cudaMalloc(&L->R_MASS_hst, 4*R_size) );
    cudaCHK( cudaMemcpy(L->R_MASS_hst, R_MASS, 4*R_size, cudaMemcpyHostToDevice) );

    cudaMalloc2E((void**)&L->R35_hst, R3, R5, R_size);
    //cudaCHK( cudaMemcpyToSymbol(R35_dev, &R35_hst, sizeof(float*)) );

    // output variables === M N Z ===
    cudaCHK( cudaMalloc(&L->Zdat_hst, L->size_hst[3]) );
    cudaCHK( cudaMalloc(&L->Zmax_hst, L->size_hst[3]) );
    //cudaCHK( cudaMemcpyToSymbol(Z_dat_dev, &Zdat_hst, sizeof(float*)) );
    cudaCHK( cudaMemcpy(L->Zdat_hst, Z_f, L->size_hst[3], cudaMemcpyHostToDevice) );
    cudaCHK( cudaMemcpy(L->Zmax_hst, Z_f, L->size_hst[3], cudaMemcpyHostToDevice) );

    cudaCHK( cudaMalloc(&L->Zout_hst, L->size_hst[3]) );
    // cudaCHK( cudaMemcpyToSymbol(Z_out_dev + L->lid, &L->Zout_hst, sizeof(float*)) ); // this'll be deleted


    cudaMalloc2E((void**)&L->MNdat_hst, NULL, NULL, L->size_hst[3]);
    //cudaCHK( cudaMemcpyToSymbol(MN_dat_dev, &MNdat_hst, sizeof(float*)) );

    cudaMalloc2E((void**)&L->MNout_hst, NULL, NULL, L->size_hst[3]);
    // cudaCHK( cudaMemcpyToSymbol(MN_out_dev + L->lid, &L->MNout_hst, sizeof(float*)) ); // this'll be deleted

    // copy data into variables
    // const variables === H, size ===
    cudaCHK( cudaMemcpy(L->H_hst, H_f, L->size_hst[3], cudaMemcpyHostToDevice) );
    //cudaCHK( cudaMemcpyToSymbol(H_dev, &H_hst, sizeof(float*)) );


    int *ptr;
    memcpy(all_size[L->lid], L->size_hst, 4*sizeof(uint32_t));
    cudaCHK( cudaMalloc(&ptr, sizeof(all_size)) );
    cudaCHK( cudaMemcpy(ptr, all_size, sizeof(all_size), cudaMemcpyHostToDevice) );
    cudaCHK( cudaMemcpyToSymbol(all_size_dev, ptr, sizeof(all_size)) );
    cudaCHK( cudaFree(ptr) );
    

    // kernel configurations
    L->DimGridMomt_MN   = dim3((L->size_hst[0]-1)/(31*(BLOCKX_MOMT>>5)) + 1, LOAD_PER_SM*(uint32_t)dev_prop.multiProcessorCount, 1);
    L->DimGridMomt      = dim3((L->size_hst[0]-1)/EXECY + 1, (L->size_hst[1]-1)/EXECX + 1, 1);
    L->DimGridMass      = dim3((L->size_hst[0]-1)/(31*(BLOCKX_MASS>>5)) + 1, LOAD_PER_SM*(uint32_t)dev_prop.multiProcessorCount, 1);
    L->DimGridOpenBD_LR = dim3((L->size_hst[0]-1)/(31*(BLOCKX_OPENBD>>5)) + 1, 1, 1);
    L->DimGridOpenBD_TB = dim3((L->size_hst[1]-1)/(31*(BLOCKX_OPENBD>>5)) + 1, 1, 1);
    L->GridMaxAmp       = (L->size_hst[2]-1)/MAXAMP_BLOCK + 1;


    #ifdef DEBUG
        // cudaError_t err;
        //
        // printf("==== R2 ele chk===\n");
        // prt_mat(R2_f, size_hst[0], size_hst[1]);
        // printf("%e\n", R2_f[ID_hst(CHKR, CHKC)]);
        //
        // printf("==== R3 ele chk===\n" );
        // prt_mat(R3_f, size_hst[0], size_hst[1]);
        // printf("%e\n", R3_f[ID_hst(CHKR, CHKC)]);
        //
        // printf("==== H ele chk===\n");
        // prt_mat(H_f, size_hst[0], size_hst[1]);
        // printf("%e\n", H_f[ID_hst(CHKR, CHKC)]);
        //
        // printf("====R5===\n");
        // prt_mat(R5_f,size_hst[0], size_hst[1]);
        //
        // printf("====Z===\n" );
        // prt_mat(Z_f,size_hst[0], size_hst[1]);
        //
        // // printf("====R6===\n");
        // // for (size_t i = CHKC; i < CHKC+CHKSI; i++) {
        // //     printf("%e\t", R6_f[i]);
        // // }
        // // printf("\n");
        // //
        // // printf("====R1====\n");
        // // prt_mat(R1_f,size_hst[0], size_hst[1]);
        // //
        // // printf("====R11====\n");
        // // prt_mat(R11_f,size_hst[0],size_hst[1]);
        //
        //
        // GPUmemtest<<< 1, 1>>>();
        // cudaDeviceSynchronize();
        // err = cudaGetLastError();
        // cudaERROR(err);
        tmpout = (float*) malloc(2*size_hst[3]);
    #endif

    free(R_MASS);
    free(R3);
    free(R5);
}

extern "C" void cuda_update_layer_(int *lid) {
    struct GPU_Layer *L = ldlayer(*lid);
    float *tmp;
    // similar to function change
    // cudaCHK( cudaMemcpyToSymbol(MN_out_dev + *lid, &MNdat_hst, sizeof(float*)) ); // need to be deleted
    //cudaCHK( cudaMemcpyToSymbol(MN_dat_dev, &MNout_hst, sizeof(float*)) );
    tmp = L->MNout_hst;
    L->MNout_hst = L->MNdat_hst;
    L->MNdat_hst = tmp;

    // cudaCHK( cudaMemcpyToSymbol(Z_out_dev + *lid, &Zdat_hst, sizeof(float*)) ); // need to be deleted
    // cudaCHK( cudaMemcpyToSymbol(Z_dat_dev, &Zout_hst, sizeof(float*)) );
    tmp = L->Zout_hst;
    L->Zout_hst = L->Zdat_hst;
    L->Zdat_hst = tmp;
}

extern "C" void gcomcot_init_gpu_(void)
{
    cudaCHK( cudaGetDeviceProperties(&dev_prop, 0) );
    printf("GPU INFORMATIONS:                      %s\n", dev_prop.name);
    printf("-->Compute Capabilities [Major.Miner]: %d.%d\n", dev_prop.major, dev_prop.minor);
    printf("-->Clock Rate:                         %d\n", dev_prop.clockRate);
    printf("-->Streaming Multi-Processor Count:    %d\n", dev_prop.multiProcessorCount);
    printf("-->Shared Memory size per SM           %d\n", dev_prop.sharedMemPerBlock);
    printf("-->Total Constant Memory size:         %d\n", dev_prop.totalConstMem );
    printf("-->Maximum Grid Size:                  %dx%dx%d\n", dev_prop.maxGridSize[0],dev_prop.maxGridSize[1],dev_prop.maxGridSize[2]);
    printf("-->Warp Size:                          %d\n", dev_prop.warpSize);
    
    //streams
    for (size_t i = 0; i < NUMSTREAM; i++) {
        cudaCHK( cudaStreamCreate(EXECstream+i) );
    }
}


extern "C" void gcomcot_init_layer_(int *layerid, int *parent, int *level,
                                
                                float *R1_f, float *R2_f,
                                float *R3_f, float *R4_f, float *R5_f, 
                                float *R6_f, float *R11_f, float *H_f, 
                                float *Z_f, int *row, int *col)
{
    struct GPU_Layer *L = ldlayer(*layerid);
    
    if (*layerid < 0) {
        printf("ERROR: invalid layerid\n");
        exit(-1);
    }
    if (*layerid >= MAX_LAYERS) {
        printf("ERROR: number of layer exceed MAX_LAYERS\n");
        exit(-1);
    }
    printf("Initializing layer id %d\n", *layerid);
    
    L->lid = *layerid;
    L->plid = *parent;
    L->lvl = *level;
    
    cuda_alloc_layer(L, R1_f, R2_f, R3_f, R4_f, R5_f,
                     R6_f, R11_f, H_f, Z_f, row, col);
    
    
}

extern "C" void cuda_shutdown_(void){
    // Free Device Memory
    // Unbind Texture Memory
}


void cudaMalloc2E(void** cu_hst, void* e1, void* e2, size_t size){
    cudaCHK( cudaMalloc(cu_hst, 2*size) );
    if (!e1) {
        cudaCHK( cudaMemset(*cu_hst, 0, size) );
    }
    else{
        cudaCHK( cudaMemcpy(*cu_hst, e1, size, cudaMemcpyHostToDevice) );
    }
    if (!e2) {
        cudaCHK( cudaMemset((void*)((char*)*cu_hst+size), 0, size) );
    }else{
        cudaCHK( cudaMemcpy((void*)((char*)*cu_hst+size), e2, size, cudaMemcpyHostToDevice) );
    }
}

#ifdef DEBUG
    __global__ void CHECK_VAR(){
        for (size_t i = 0; i < size_dev[0]; i++) {
            printf("%e\t",Z_out_dev[i]);
        }
        printf("\n" );
    }
    extern "C" void cmpvar_(const float *var_f, int *Case){
        if (*Case == 0) {
            printf("copying data\n");
            for (size_t i = 0; i < size_hst[2]; i++) {
                tmpout[i] = var_f[i];
            }
            return;
        }else{
            printf("comparing data\n" );
            for (size_t i = 0; i < size_hst[2]; i++) {
                if (abs(tmpout[i] - var_f[i]) > ERROR) {
                    printf("VAR[%d,%d] VAR_cu:%e VAR_f:%e %e\n", i%size_hst[0], i/size_hst[0] , tmpout[i], var_f[i], tmpout[i] - var_f[i]);
                }
            }
        }
    }
    // __global__ void GPUmemtest() {
    //     printf("%d %d %d %d\n",size_dev[0], size_dev[1], size_dev[2], size_dev[3]);
    //     printf("xxxxx GPU H ele chk xxxxx\n" );
    //     prt_mat(H_dev, size_dev[0], size_dev[1]);
    //     printf("%f\n", H_dev[ID(CHKR, CHKC)]);
    //
    //     printf("xxxxx GPU R2 ele chk xxxxx\n" );
    //     prt_mat(R24_dev, size_dev[0], size_dev[1]);
    //     printf("%f\n", R24_dev[ID(CHKR, CHKC)]);
    //
    //     printf("xxxxx GPU R3 ele chk xxxxx\n");
    //     for (size_t i = CHKC; i < CHKC+CHKSI; i++) {
    //         printf("%8.6e\t", R35_dev[i]);
    //     }
    // }


    __host__  __device__ void prt_mat(float *mat, size_t row, size_t col){
        printf("Matrix[%d:%d][%d:%d]\n", CHKR,CHKR+CHKSI, CHKC,CHKC+CHKSI);
        for (size_t i = CHKR; i < CHKR+CHKSI; i++) {
            for (size_t j = CHKC; j < CHKC+CHKSI; j++) {
                printf("%8.6e\t",  mat[j*row + i]);
            }
            printf("\n");
        }
    }
#endif
