#ifndef SM_CONFG
    #define SM_CONFG
    #define LOAD_PER_SM 2
    #define NUMSTREAM 4
    extern cudaStream_t EXECstream[NUMSTREAM];
#endif

#ifndef MOMT_KERNEL_CONFG
    #define MOMT_KERNEL_CONFG
    #define BLOCKX 16 // ==> along column axis
    #define EXECX  15 // BLOCKX-1
    #define BLOCKY 16 // ==> along row axis
    #define EXECY  15// BLOCKY-1

     extern dim3 DimBlockMomt;
     extern dim3 DimGridMomt;

     #define BLOCKX_MOMT 512
     extern dim3 DimBlockMomt_MN;
     extern dim3 DimGridMomt_MN;
#endif

#ifndef MASS_KERNEL_CONFG
    #define MASS_KERNEL_CONFG
    #define BLOCKX_MASS 64
     extern dim3 DimBlockMass;
     extern dim3 DimGridMass;
#endif

#ifndef OPENBD_KERNEL_CONFG
    #define OPENBD_KERNEL_CONFG
    #define BLOCKX_OPENBD 64
    extern dim3 DimBlockOpenBD;
    extern dim3 DimGridOpenBD_LR;
    extern dim3 DimGridOpenBD_TB;
#endif

#ifndef MAXAMP_KERNEL_CONFG
    #define MAXAMP_KERNEL_CONFG
    #define MAXAMP_BLOCK 512
    extern size_t GridMaxAmp;
    #define BLOCKX_MAXAMP 32
    #define BLOCKY_MAXAMP 32
    extern dim3 DimBlockMaxAmp;
    extern dim3 DimGridMaxAmp;
#endif
