// data structure for efficiently load and store data.
// It contains a 4x4 half tile of unsigned chars to keep the values of booleans.
// Since the kernels update data in chessboard manner the 4x4 tiles are divided in
// black and white elements, this means that the data structure only contains either black
// for white elements. In order to have all elements of the entire 4x4 tile you need two structures.
#include "uintOfBits.h"
#include <cuda_runtime.h>
#include <stdint.h>
#include <stdio.h>

struct CubeRedElements { // waisting 4 bits per cube
    union {
        u_int8_t cube : 4;
        struct {
            u_int8_t el_tlf : 1; // top-left-front
            u_int8_t el_trb : 1; // top-right-back
            u_int8_t el_blb : 1; // bottom-left-back
            u_int8_t el_brf : 1; // bottom-right-front
        };
    };
};

template <typename CubeColourElements> struct Tile {
    union {
        uint_of_bits<8 * sizeof(CubeColourElements)> tile;
        struct {
            CubeColourElements c_tlf; // cube top-left-front
            CubeColourElements c_tlb; // cube top-left-back
            CubeColourElements c_trf; // cube top-right-front
            CubeColourElements c_trb; // cube top-right-back

            CubeColourElements c_blf; // cube bottom-left-front
            CubeColourElements c_blb; // cube bottom-left-back
            CubeColourElements c_brf; // cube bottom-right-front
            CubeColourElements c_brb; // cube bottom-right-back
        };
    };
};

__global__ void efficientKernel(Tile<CubeRedElements> *matrix, int Ndiv4) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < Ndiv4) {
        Tile tile = matrix[idx];
        uint a = sizeof(tile.tile);
        u_int32_t c_tlf_el_brb = (u_int32_t)(tile.c_tlf.el_blb) + (u_int32_t)tile.c_tlf.el_brf +
                                 (u_int32_t)tile.c_tlf.el_trb + (u_int32_t)tile.c_tlb.el_brf +
                                 (u_int32_t)tile.c_blf.el_trb + (u_int32_t)tile.c_trf.el_blb;
        volatile double x = 0.11 * c_tlf_el_brb;
    }
}

int main(void) {

    constexpr int N = 1024; // Example size
    constexpr int numTiles = N / 4;

    Tile<CubeRedElements> *d_matrix;
    size_t size = numTiles * sizeof(Tile<CubeRedElements>);
    cudaMalloc((void **)&d_matrix, size);
    cudaMemset(d_matrix, 0, size); // Initialize the device memory
    int threadsPerBlock = numTiles;
    int blocksPerGrid = (numTiles + threadsPerBlock - 1) / threadsPerBlock;

    efficientKernel<<<1, threadsPerBlock>>>(d_matrix, numTiles);
    cudaDeviceSynchronize();
    cudaFree(d_matrix);
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(err));
        return -1;
    }
    printf("Kernel executed successfully.\n");

    return 0;
}
