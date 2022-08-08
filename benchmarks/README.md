# Benchmarks
Current file is a similar to our electronic structure calculation script.
For better performance do not use `total number of cores` that your computer provides. Instead, try to use `total number of cores - 1` as general. Use `time` command to measure the time passed as (prefered version of running GPAW is):


       gpaw -P7 python simple_benchmark_2021.py

or (this will take a little bit longer)

       time mpiexec -n 7 gpaw python simple_benchmark_2021.py
       
or if your CPU supports threads (sometimes this may take even longer)

       time mpiexec --use-hwthread-cpus -n 15 gpaw python simple_benchmark_2021.py

These commands will result something like:
```
       Step     Time          Energy         fmax
*Force-consistent energies used in optimization.
LBFGS:    0 14:20:09     -141.676163*       1.6416
LBFGS:    1 14:20:59     -141.790250*       0.9145
LBFGS:    2 14:21:37     -141.860481*       0.7577
LBFGS:    3 14:22:37     -142.187974*       0.7607
LBFGS:    4 14:23:13     -142.253337*       0.4721
LBFGS:    5 14:23:49     -142.283340*       0.1869
LBFGS:    6 14:24:27     -142.286176*       0.1619
LBFGS:    7 14:25:00     -142.289792*       0.1506
LBFGS:    8 14:25:35     -142.293183*       0.1230
LBFGS:    9 14:26:02     -142.294811*       0.0715
LBFGS:   10 14:26:33     -142.295239*       0.0274

real    9m0.719s
user    62m4.101s
sys     0m40.817s
```
Here, `real    9m0.0719s` is the benchmark time.

## Some Benchmark Times

### Computers
| Computer  | CPU                      | Cores(Threads) | CPU Speed | Memory | Hdd            | Sysbench Events(1) | Memory Bandwidth | HDD Speed |
| --------- | ------------------------ | -------------- | --------- | ------ | -------------- | ------------------ | ---------------- | --------- |
| 1         | 2x Intel Xeon E5-2430 v2 | 12(24)         | 2.5 GHz   | 16Gb   | 300Gb + 1000Gb | 15525.13           | 5 GB/s           | 259 MB/s  |
| 2         | Intel Core i7-8550u      | 4(8)           | 1.8 GHz   | 8Gb    | 512SSD         | -                  | 12 GB/s          | 1GB/s     |
| 3         | AMD Ryzen 5 4500u        | 6(6)           | 2.4 GHz   | 8Gb    | 256SSD         | -                  | 9 GB/s           | 0.8GB/s   |
| 4         | Intel Core i7-8700       | 6(12)          | 3.2 GHz   | 8Gb    | 1000 Gb        | 13302.23           | 10 GB/s          | 195 MB/s  |
| 5         | AMD Ryzen 7 5700u        | 8(16)          | 1.8 GHz   | 16Gb   | 512SSD         | 17201.45           | 6 GB/s           | 2GB/s     |
| 6         | AMD Ryzen 7 4700G        | 8(16)          | 3.6 GHz   | 16Gb   | 512SSD         | 18466.64           | 10 GB/s          | 1 GB/s    |
| 7         | Intel Xeon W3540         | 4(8)           | 2.9 GHz   | 9Gb    | -              | -                  | -                | -         |
| 8         | Intel i5-8250u           | 4(8)           | 1.6 GHz   | 8 Gb   | -              | -                  | -                | -         |

### Benchmarks
| Computer  | GPAW Version  | System                  | Used Core | Command                       | Benchmark File             | Time Elapsed |
| --------- | ------------- | ----------------------- | --------- | ----------------------------- | -------------------------- | ------------ |
| 1         | 21.6.0        | W10Pro - WSL1 - Ub20.04 | 12        | gpaw -P                       | simple_benchmark_2021.py   | 5m41s        |
| 1         | 21.6.0        | W10Pro - WSL2 - Ub20.04 | 12        | gpaw -P                       | simple_benchmark_2021.py   | 5m59s        |
| 1         | 21.6.0        | W10Pro - WSL1 - Ub20.04 | 23        | mpirun --use-hwthread-cpus -n | simple_benchmark_2021.py   | 5m15s        |
| 1         | 21.6.0        | W10Pro - WSL2 - Ub20.04 | 23        | mpirun --use-hwthread-cpus -n | simple_benchmark_2021.py   | 5m37s        |
| 2         | 21.6.0        | W10Pro - WSL1 - Ub20.04 | 7         | mpirun --use-hwthread-cpus -n | simple_benchmark_2021.py   | 7m24s        |
| 2         | 21.6.0        | W10Pro - WSL2 - Ub20.04 | 7         | mpirun --use-hwthread-cpus -n | simple_benchmark_2021.py   | 9m00s        |
| 3         | 21.6.0        | W10 - WSL1 - Ub20.04    | 5         | gpaw -P                       | simple_benchmark_2021.py   | 4m26s        |
| 4         | 21.6.0        | W10Pro - WSL1 - Ub20.04 | 11        | mpirun --use-hwthread-cpus -n | simple_benchmark_2021.py   | 5m22s        |
| 5         | 21.6.0        | Ubuntu 20.04            | 8         | gpaw -P                       | simple_benchmark_2021.py   | 4m24s        |
| 5         | 21.6.0        | Ubuntu 20.04            | 8         | mpirun -n                     | simple_benchmark_2021.py   | 5m15s        |
| 5         | 21.6.0        | Ubuntu 20.04            | 15        | mpirun --use-hwthread-cpus -n | simple_benchmark_2021.py   | 8m19s        |
| 6         | 21.6.0        | W10Pro - WSL1 - Ub20.04 | 8         | gpaw -P                       | simple_benchmark_2021.py   | 3m50s        |
| 7         | 21.6.0        | W10Pro - WSL1 - Ub20.04 | 8         | mpirun -n                     | simple_benchmark_2021.py   | 11m25s       |
| 8         | 22.1.0        | W10Pro - WSL1 - Ub20.04 | 2         | mpirun -n                     | simple_benchmark_2021.py   | 15m28s       |

(1) : sysbench for linux is used for CPU benchmark (UBUNTU INSTALLATION: sudo apt-get install sysbench | USAGE: sysbench --test=cpu --threads=PUT-THREAD-NUMBER-HERE run).
