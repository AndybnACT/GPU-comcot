#ifndef GPUCONF_H
#define GPUCONF_H

#define SM_CONFG
#define LOAD_PER_SM 2
#define NUMSTREAM 4
extern cudaStream_t EXECstream[NUMSTREAM];

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

#define BLOCKX_MASS 64
extern dim3 DimBlockMass;
extern dim3 DimGridMass;

#define BLOCKX_JNZ 64
#define DEBUG_DEPTH 32
#define LOOPDEPTH DEBUG_DEPTH

#define BLOCKX_OPENBD 64
extern dim3 DimBlockOpenBD;
extern dim3 DimGridOpenBD_LR;
extern dim3 DimGridOpenBD_TB;

#define MAXAMP_BLOCK 512
extern size_t GridMaxAmp;

#define BLOCKX_JNQ 128

#endif /* GPUCONF_H */
