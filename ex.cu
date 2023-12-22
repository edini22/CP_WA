#include <iostream>
#include <cstdlib>
#include <cuda.h>

// Função do kernel CUDA para somar dois vetores
__global__ void addVectors(int *a, int *b, int *c, int size) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid < size) {
        c[tid] = a[tid] + b[tid];
    }
}

// Função que realiza a chamada da função CUDA
void performCUDAOperation(int *d_a, int *d_b, int *d_c, int size) {
    // Configurar a grade de threads e os blocos
    int threadsPerBlock = 256;
    int blocksPerGrid = (size + threadsPerBlock - 1) / threadsPerBlock;

    // Chamar o kernel CUDA para somar os vetores
    addVectors<<<blocksPerGrid, threadsPerBlock>>>(d_a, d_b, d_c, size);

    // Sincronizar o dispositivo para garantir que o kernel seja concluído
    cudaDeviceSynchronize();

    // Copiar o resultado de volta para o host (aqui você pode realizar operações adicionais se necessário)
    // Exemplo: cudaMemcpy(c, d_c, size * sizeof(int), cudaMemcpyDeviceToHost);
}

int main() {
    const int size = 10; // Tamanho dos vetores
    const int iterations = 5; // Número de iterações do loop
    int a[size], b[size], c[size]; // Vetores de entrada e saída no host
    int *d_a, *d_b, *d_c; // Vetores no device (GPU)

    // Inicialização dos vetores no host
    for (int i = 0; i < size; ++i) {
        a[i] = i;
        b[i] = 2 * i;
    }

    // Alocação de memória no device
    cudaMalloc((void **)&d_a, size * sizeof(int));
    cudaMalloc((void **)&d_b, size * sizeof(int));
    cudaMalloc((void **)&d_c, size * sizeof(int));

    // Copiar dados do host para o device
    cudaMemcpy(d_a, a, size * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, b, size * sizeof(int), cudaMemcpyHostToDevice);

    // Loop para realizar várias iterações
    for (int iter = 0; iter < iterations; ++iter) {
        // Chamar a função que executa a operação CUDA
        performCUDAOperation(d_a, d_b, d_c, size);

        // Copiar os resultados de volta para o host
        cudaMemcpy(c, d_c, size * sizeof(int), cudaMemcpyDeviceToHost);

        // Modificar as variáveis no host com base nos resultados do kernel
        for (int i = 0; i < size; ++i) {
            a[i] += 1;  // Incrementar 1 em a com base nos resultados do kernel
        }

        // Copiar as variáveis modificadas de volta para o device (se necessário)
        cudaMemcpy(d_a, a, size * sizeof(int), cudaMemcpyHostToDevice);
    }

    // Exibir o resultado final (aqui você pode realizar operações adicionais se necessário)
    for (int i = 0; i < size; i++)
    {
        std::cout << a[i] << " ";
    }
    std::cout << "\n";
    for (int i = 0; i < size; i++)
    {
        std::cout << b[i] << " ";
    }
    std::cout << "\n";
    std::cout << "Resultado da soma final:\n";
    for (int i = 0; i < size; ++i) {
        std::cout << c[i] << " ";
    }
    std::cout << "\n";
    // Liberar memória no device
    cudaFree(d_a);
    cudaFree(d_b);
    cudaFree(d_c);

    return 0;
}
