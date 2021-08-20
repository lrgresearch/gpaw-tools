# Benchmarks
Current file is a similar to our electronic structure calculation script.
For better performance do not use `total number of cores` that your computer provides. Instead, try to use `total number of cores - 1` as general. Use `time` command to measure the time passed as (prefered version of running GPAW is):

```
gpaw -P7 python GPAWSimpleBenchmark2021.py
```
or (this will take a little bit longer)
```
time mpiexec -n 7 gpaw python GPAWSimpleBenchmark2021.py
```
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
| Computer  | CPU                      | Cores | CPU Speed | Memory | Hdd            | CPU GFlops | Memory Bandwidth | HDD Speed |
| --------- | ------------------------ | ----- | --------- | ------ | -------------- | ---------- | ---------------- | --------- |
| 1         | 2x Intel Xeon E5-2430 v2 | 24    | 2.5 GHz   | 16Gb   | 300Gb + 1000Gb | 357        | 5 GB/s           | 259 MB/s  |

### Benchmarks
| Computer  | GPAW Version  | System                  | Used Core | Command                       | Benchmark File             | Time Elapsed |
| --------- | ------------- | ----------------------- | --------- | ----------------------------- | -------------------------- | ------------ |
| 1         | 21.6.0        | W10Pro - WSL1 - Ub20.04 | 12        | gpaw -P                       | GPAWSimpleBenchmark2021.py | 5m41s        |
| 1         | 21.6.0        | W10Pro - WSL2 - Ub20.04 | 12        | gpaw -P                       | GPAWSimpleBenchmark2021.py | 5m59s        |
| 1         | 21.6.0        | W10Pro - WSL1 - Ub20.04 | 23        | mpirun --use-hwthread-cpus -n | GPAWSimpleBenchmark2021.py | 5m15s        |
| 1         | 21.6.0        | W10Pro - WSL1 - Ub20.04 | 23        | mpirun --use-hwthread-cpus -n | GPAWSimpleBenchmark2021.py | 5m37s        |

