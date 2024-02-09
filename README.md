# Optimizing Molecular Dynamics Simulation for Argon Gas

## Authors
- [Eduardo Figueiredo]
- [Gonçalo Senra]

## Contacts
- [Eduardo Figueiredo]
    - Email: pg52679@uminho.pt

- [Gonçalo Senra]
    - Email: pg52683@uminho.pt

# Abstract:
This article explores techniques to enhance the performance of a computer program simulating the behavior of argon gas. Specialized optimization methods and parallel computing paradigms are employed to improve the program's efficiency and speed. The study includes the analysis of sequential code optimization, parallelization using OpenMP, and optimization with CUDA for GPU acceleration. Results and insights gained from various tests and implementations are presented to provide a comprehensive understanding of the optimization process.

# Introduction:
The focus of this study is to analyze and enhance a computer program responsible for simulating the behavior of argon gas molecules. Initially, efforts were directed towards optimizing the sequential execution of the program by restructuring code and improving algorithm efficiency. Subsequently, the program was parallelized using OpenMP to leverage multicore processors. Finally, CUDA was employed to achieve further optimization by harnessing the parallel processing capabilities of GPUs. The methodologies employed in each phase aimed to significantly enhance the program's performance and efficiency.

# Sequential Code Optimization:
The optimization process involved several strategies aimed at improving the efficiency of the molecular dynamics simulation code. These strategies included optimizing potential calculation, accelerating computations, eliminating unnecessary loops, and optimizing power operations. By refining these aspects of the code, significant improvements in execution time and performance were achieved.

# Parallelization with OpenMP:
Parallelization using OpenMP involved identifying hotspots in the code and implementing parallel execution using directives such as #pragma omp parallel and #pragma omp for. By distributing computational tasks across multiple threads, the program's performance was enhanced, leading to faster execution times. Various optimization techniques and considerations were employed to ensure efficient parallel execution and avoid data race conditions.

# Optimization with CUDA:
Utilizing CUDA for GPU acceleration enabled further optimization of the simulation code. By parallelizing critical sections of the code and leveraging the parallel processing power of GPUs, significant improvements in execution time were achieved. The implementation involved creating GPU kernels, optimizing memory access, and addressing numerical precision considerations to ensure accurate results.

# Testing Methodology and Performance Evaluation:
Extensive testing was conducted to evaluate the performance of the optimized code under different conditions and input sizes. Results were analyzed to identify the impact of parallelization techniques and hardware architecture on execution time and efficiency. Insights gained from testing were used to refine implementation strategies and optimize code parameters for improved performance.

# Conclusion:
The optimization efforts undertaken in this study resulted in significant improvements in the performance and efficiency of the molecular dynamics simulation program for argon gas. By employing sequential code optimization, parallelization with OpenMP, and optimization with CUDA, the program's execution time was greatly reduced, enhancing its practical utility and scalability. The insights gained from this study provide valuable guidance for future optimization efforts in computational science and parallel computing applications.

To view the complete report of this work, it can be found [report](/code_wa3/Report_CP_WA3.pdf).
