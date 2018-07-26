#include "GPUHeader.h"
#include "GPUConfig.h"

// global variables-----------------------------
float *Zout_hst/*, *MNout_hst*/;
float /**MNdat_hst,*/ *Zdat_hst;
float *R24_hst, *R35_hst, *H_hst;
float *R_MASS_hst;
float *Zmax_hst;
// __device__ float *R35_dev;
// __device__ float *R24_dev, *H_dev;
// __device__ float *Z_dat_dev, *MN_dat_dev;
__device__ float /**MN_out_dev,*/ *Z_out_dev;
__constant__ __device__ uint32_t size_dev[4];
// texture<float, cudaTextureType2D, cudaReadModeElementType> ZtexRef;
// texture<float, cudaTextureType2D, cudaReadModeElementType> MNtexRef;

float2 *MNcontainer;
texture <float2, cudaTextureType2D, cudaReadModeElementType> MNtext;
size_t MNpitch;
float2 *MNdat_pitchedMEM_hst;
float2 *MNout_pitchedMEM_hst;

float2 *ZHcontainer;
texture <float2, cudaTextureType2D, cudaReadModeElementType> ZHtext;
size_t ZHpitch;
float2 *ZHdat_pitchedMEM_hst;
float2 *ZHout_pitchedMEM_hst;

cudaChannelFormatDesc descflt2;

__device__ float2 *MNout_pitchedMEM_dev;
__constant__ __device__ size_t MNpitch_dev;
__device__ float2 *ZHout_pitchedMEM_dev;
__constant__ __device__ size_t ZHpitch_dev;
// __constant__ size_t textOffset;

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

uint32_t size_hst[4];
cudaDeviceProp dev_prop;
#ifdef DEBUG
    float *tmpout;
#endif
//----------------------------------------------
extern "C" void cuda_boot_(float*,float*,float*,float*,float*,float*,float*,float*,float*,int*,int*);
extern "C" void cuda_update_(void);
void cudaMalloc2E(void**, void*, void*, size_t);
void cudaMallocMNZH(size_t, size_t, float*, float*);
#ifdef DEBUG
    __host__ __device__ void prt_mat(float*, size_t, size_t);
    __global__ void GPUmemtest(void);
    __global__ void CHECK_VAR(void);
    extern "C" void cmpvar_(const float *, int *);
#endif


extern "C" void cuda_boot_(float *R1_f, float *R2_f, float *R3_f, float *R4_f, float *R5_f, \
                           float *R6_f, float *R11_f, float *H_f, float *Z_f, int *nx, int *ny){

    //float *R24_hst, *R35_hst, *H_hst;
    float *R3, *R5;
    float *R_MASS;



    size_hst[0] = *nx;
    size_hst[1] = *ny;
    size_hst[2] = (*nx)*(*ny);
    size_hst[3] = (*nx)*(*ny)*sizeof(float);
    cudaCHK( cudaGetDeviceProperties(&dev_prop, 0) );
    printf("GPU INFORMATIONS:                      %s\n", dev_prop.name);
    printf("-->Compute Capabilities [Major.Miner]: %d.%d\n", dev_prop.major, dev_prop.minor);
    printf("-->Clock Rate:                         %d\n", dev_prop.clockRate);
    printf("-->Streaming Multi-Processor Count:    %d\n", dev_prop.multiProcessorCount);
    printf("-->Shared Memory size per SM           %d\n", dev_prop.sharedMemPerBlock);
    printf("-->Total Constant Memory size:         %d\n", dev_prop.totalConstMem );
    printf("-->Maximum Grid Size:                  %dx%dx%d\n", dev_prop.maxGridSize[0],dev_prop.maxGridSize[1],dev_prop.maxGridSize[2]);
    printf("-->Warp Size:                          %d\n", dev_prop.warpSize);
    // allocate variables
    // const === R2 R4 H ===
    cudaCHK( cudaMalloc(&H_hst, size_hst[3]) );

    cudaMalloc2E((void**)&R24_hst, R2_f, R4_f, size_hst[3]);
    //cudaCHK( cudaMemcpyToSymbol(R24_dev, &R24_hst, sizeof(float*)) );
    // const variables R3 R5, R1 R6 R11
    size_t R_size = sizeof(float)*size_hst[1];
    //R1 = (float*) malloc(R_size);
    R3 = (float*) malloc(R_size);
    R5 = (float*) malloc(R_size);
    //R11 =(float*) malloc(R_size);
    R_MASS = (float*) malloc(4*R_size);

    for (size_t i = 0; i < size_hst[1]; i++) {
        //R1[i] = R1_f[i*size_hst[0]];
        R3[i] = R3_f[i*size_hst[0]];
        R5[i] = R5_f[i*size_hst[0]];
        //R11[i]= R11_f[i*size_hst[0]];
        R_MASS[4*i]   = R1_f[i*size_hst[0]];
        R_MASS[4*i+1] = R6_f[i];
        R_MASS[4*i+2] = R11_f[i*size_hst[0]];
    }


    // cudaCHK( cudaMalloc(&R1_hst, R_size) );
    // cudaCHK( cudaMemcpy(R1_hst, R1, R_size, cudaMemcpyHostToDevice) );

    // cudaCHK( cudaMalloc(&R6_hst, R_size) );
    // cudaCHK( cudaMemcpy(R6_hst, R6_f, R_size, cudaMemcpyHostToDevice) );

    // cudaCHK( cudaMalloc(&R11_hst, R_size) );
    // cudaCHK( cudaMemcpy(R11_hst, R11, R_size, cudaMemcpyHostToDevice) );
    cudaCHK( cudaMalloc(&R_MASS_hst, 4*R_size) );
    cudaCHK( cudaMemcpy(R_MASS_hst, R_MASS, 4*R_size, cudaMemcpyHostToDevice) );

    cudaMalloc2E((void**)&R35_hst, R3, R5, R_size);
    //cudaCHK( cudaMemcpyToSymbol(R35_dev, &R35_hst, sizeof(float*)) );

    // output variables === M N Z ===
    cudaCHK( cudaMalloc(&Zdat_hst, size_hst[3]) );
    cudaCHK( cudaMalloc(&Zmax_hst, size_hst[3]) );
    //cudaCHK( cudaMemcpyToSymbol(Z_dat_dev, &Zdat_hst, sizeof(float*)) );
    cudaCHK( cudaMemcpy(Zdat_hst, Z_f, size_hst[3], cudaMemcpyHostToDevice) );
    cudaCHK( cudaMemcpy(Zmax_hst, Z_f, size_hst[3], cudaMemcpyHostToDevice) );

    cudaCHK( cudaMalloc(&Zout_hst, size_hst[3]) );
    cudaCHK( cudaMemcpyToSymbol(Z_out_dev, &Zout_hst, sizeof(float*)) );


    // cudaMalloc2E((void**)&MNdat_hst, NULL, NULL, size_hst[3]);
    //cudaCHK( cudaMemcpyToSymbol(MN_dat_dev, &MNdat_hst, sizeof(float*)) );

    // cudaMalloc2E((void**)&MNout_hst, NULL, NULL, size_hst[3]);
    // cudaCHK( cudaMemcpyToSymbol(MN_out_dev, &MNout_hst, sizeof(float*)) );

    // copy data into variables
    // const variables === H, size ===
    cudaCHK( cudaMemcpy(H_hst, H_f, size_hst[3], cudaMemcpyHostToDevice) );
    //cudaCHK( cudaMemcpyToSymbol(H_dev, &H_hst, sizeof(float*)) );

    cudaCHK( cudaMemcpyToSymbol(size_dev, &size_hst, sizeof(size_hst)) );

    //Texture Memories
    // cudaChannelFormatDesc desc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
    // cudaCHK( cudaBindTexture2D(NULL, MNtext, M) );
    // MN texture
    cudaMallocMNZH(*ny, *nx, Z_f, H_f);
    descflt2 = cudaCreateChannelDesc(32, 32, 0, 0, cudaChannelFormatKindFloat);
    MNtext.filterMode = cudaFilterModePoint;
    MNtext.normalized = 0;
    MNtext.addressMode[0] = cudaAddressModeClamp;
    MNtext.addressMode[1] = cudaAddressModeClamp;

    ZHtext.filterMode = cudaFilterModePoint;
    ZHtext.normalized = 0;
    ZHtext.addressMode[0] = cudaAddressModeClamp;
    ZHtext.addressMode[1] = cudaAddressModeClamp;
    //                                                                       |-> width=col  |-> height=row
    cudaCHK( cudaBindTexture2D(NULL, MNtext, MNdat_pitchedMEM_hst, descflt2, size_hst[0], size_hst[1], MNpitch) );
    cudaCHK( cudaBindTexture2D(NULL, ZHtext, ZHdat_pitchedMEM_hst, descflt2, size_hst[0], size_hst[1], MNpitch) );


    // kernel configurations
    DimGridMomt_MN   = dim3((size_hst[0]-1)/(31*(BLOCKX_MOMT>>5)) + 1, LOAD_PER_SM*(uint32_t)dev_prop.multiProcessorCount, 1);
    DimGridMomt      = dim3((size_hst[0]-1)/EXECY + 1, (size_hst[1]-1)/EXECX + 1, 1);
    DimGridMass      = dim3((size_hst[0]-1)/(31*(BLOCKX_MASS>>5)) + 1, LOAD_PER_SM*(uint32_t)dev_prop.multiProcessorCount, 1);
    DimGridOpenBD_LR = dim3((size_hst[0]-1)/(31*(BLOCKX_OPENBD>>5)) + 1, 1, 1);
    DimGridOpenBD_TB = dim3((size_hst[1]-1)/(31*(BLOCKX_OPENBD>>5)) + 1, 1, 1);
    GridMaxAmp       = (size_hst[2]-1)/MAXAMP_BLOCK + 1;

    //streams
    for (size_t i = 0; i < NUMSTREAM; i++) {
        cudaCHK( cudaStreamCreate(EXECstream+i) );
    }

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

extern "C" void cuda_update_(void) {// similar to function change
    float *tmp;
    float2 *tmpf2;
    // unbind the MN texture
    cudaCHK( cudaUnbindTexture(MNtext) );
    cudaCHK( cudaUnbindTexture(ZHtext) );
    // FUTURE developments should NOTE that
    // ZHtext here has already been binded to ZHout_pitchedMEM_hst by momt_launch_
    // since momts_kernel_tex required the updated Z.


    // cudaCHK( cudaMemcpyToSymbol(MN_out_dev, &MNdat_hst, sizeof(float*)) );
    // //cudaCHK( cudaMemcpyToSymbol(MN_dat_dev, &MNout_hst, sizeof(float*)) );
    // tmp = MNout_hst;
    // MNout_hst = MNdat_hst;
    // MNdat_hst = tmp;

    cudaCHK( cudaMemcpyToSymbol(Z_out_dev, &Zdat_hst, sizeof(float*)) );
    // cudaCHK( cudaMemcpyToSymbol(Z_dat_dev, &Zout_hst, sizeof(float*)) );
    tmp = Zout_hst;
    Zout_hst = Zdat_hst;
    Zdat_hst = tmp;

    cudaCHK( cudaMemcpyToSymbol(MNout_pitchedMEM_dev, &MNdat_pitchedMEM_hst, sizeof(float*)) );
    tmpf2 = MNout_pitchedMEM_hst;
    MNout_pitchedMEM_hst = MNdat_pitchedMEM_hst;
    MNdat_pitchedMEM_hst = tmpf2;

    cudaCHK( cudaMemcpyToSymbol(ZHout_pitchedMEM_dev, &ZHdat_pitchedMEM_hst, sizeof(float*)) );
    tmpf2 = ZHout_pitchedMEM_hst;
    ZHout_pitchedMEM_hst = ZHdat_pitchedMEM_hst;
    ZHdat_pitchedMEM_hst = tmpf2;

    // rebind the texture
    //                                                                        |-> width=col  |-> height=row
    cudaCHK( cudaBindTexture2D(NULL, MNtext, MNdat_pitchedMEM_hst, descflt2, size_hst[0], size_hst[1], MNpitch) );
    cudaCHK( cudaBindTexture2D(NULL, ZHtext, ZHdat_pitchedMEM_hst, descflt2, size_hst[0], size_hst[1], MNpitch) );

}

extern "C" void cuda_shutdown_(void){
    // Free Device Memory
    // Unbind Texture Memory
}


void cudaMallocMNZH(size_t ny, size_t nx, float* Z_f, float* H_f){
    size_t tmppitch;
    MNcontainer = (float2 *) malloc(nx*ny*sizeof(float2));
    ZHcontainer = (float2 *) malloc(nx*ny*sizeof(float2));
    for (size_t i = 0; i < ny; i++) {// assign initial value to host with float2 format
        for (size_t j = 0; j < nx; j++) {
            // MNcontainer[i*nx + j].x = 0.0f; //use memset instead
            // MNcontainer[i*nx + j].y = 0.0f;
            ZHcontainer[i*nx + j].x = Z_f[i*nx + j];
            ZHcontainer[i*nx + j].y = H_f[i*nx + j];
        }
    }

    //device memory--dat
    cudaCHK( cudaMallocPitch(&MNdat_pitchedMEM_hst, &MNpitch, nx*sizeof(float2), ny) );
    cudaCHK( cudaMallocPitch(&ZHdat_pitchedMEM_hst, &tmppitch, nx*sizeof(float2), ny) );
    if (tmppitch != MNpitch) { // just a check
        printf("Error @%d, %s\n",__LINE__, __FILE__ ); // not reach
        exit(1);
    }
    //device memory--out
    cudaCHK( cudaMallocPitch(&MNout_pitchedMEM_hst, &tmppitch, nx*sizeof(float2), ny) );
    cudaCHK( cudaMallocPitch(&ZHout_pitchedMEM_hst, &tmppitch, nx*sizeof(float2), ny) );

    //device symbol
    cudaCHK( cudaMemcpyToSymbol(MNout_pitchedMEM_dev, &MNout_pitchedMEM_hst, sizeof(float2*)) );
    cudaCHK( cudaMemcpyToSymbol(ZHout_pitchedMEM_dev, &ZHout_pitchedMEM_hst, sizeof(float2*)) );

    //pitch symbol
    cudaCHK( cudaMemcpyToSymbol(MNpitch_dev, &MNpitch, sizeof(size_t)) );

    // copy the initial values into device
    cudaCHK( cudaMemset2D(MNdat_pitchedMEM_hst, MNpitch, 0, nx*sizeof(float2), ny) );
    // cudaCHK( cudaMemcpy2D(MNdat_pitchedMEM_hst, MNpitch, MNcontainer, nx*sizeof(float2)
    //                                             , nx*sizeof(float2), ny, cudaMemcpyHostToDevice) );

    cudaCHK( cudaMemcpy2D(ZHdat_pitchedMEM_hst, MNpitch, ZHcontainer, nx*sizeof(float2)//since cuda_update_ switches out/dat pointers
                                                , nx*sizeof(float2), ny, cudaMemcpyHostToDevice) ); // , H must exist in both memory
    cudaCHK( cudaMemcpy2D(ZHout_pitchedMEM_hst, MNpitch, ZHcontainer, nx*sizeof(float2)//since cuda_update_ switches out/dat pointers
                                                , nx*sizeof(float2), ny, cudaMemcpyHostToDevice) ); // , H must exist in both memory


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
