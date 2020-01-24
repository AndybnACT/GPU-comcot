#ifndef GPUOpen_BD_H
    #define GPUOpen_BD_H

    typedef enum bd_side{
        LEFT,
        RIGHT,
        TOP,
        BOTTOM,
    } bdside;


    extern "C" void openbd_launch_(float*);
    __global__ void openbd_kernel(struct GPU_Layer, bdside);


#endif
