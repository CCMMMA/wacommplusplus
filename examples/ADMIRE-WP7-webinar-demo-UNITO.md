## Start building the code

1) Clone the repository
```bash
git clone https://github.com/CCMMMA/wacommplusplus.git
```
2) Enter the project directory
```bash
cd wacomplusplus
```
3) Create the build directory
```bash
mkdir build
```
4) Enter the build directory
```bash
cd build
```
5) Invoke cmake using the following options:
- USE_MPI - Activate the distributed memory parallelization.
- USE_OMP - Activate the shared memory parallelization.
- USE_EMPI - Activate computational malleability.
- USE_CUDA - Activate the CUDA acceleration.
- DEBUG - Add logging printouts (do not use for production or evaluation.)

(WaComM++ uses cmake version 3. In some Linux systems, cmake is version 3 by default. In other systems version 3
must be specified using the cmake3 command. )

- Building on the UNITO server requires setting the environment:
```bash
spack load cmake@3.26.3%gcc@=9.4.0 arch=linux-ubuntu20.04-broadwell
spack load openmpi@4.1.5%gcc@=9.4.0 arch=linux-ubuntu20.04-zen
```

- We are going to compile with OpenMP and MPI support but without extra debug logging:
```bash
cmake -DUSE_OMP=ON -DUSE_MPI=ON -DUSE_EMPI=OFF -DUSE_CUDA=OFF -DDEBUG=OFF ..
```

6) Run, make, and wait

```bash
make
```

# Testing
Download the data.

1) Open a console application (I.e. Terminal)
2) Change the current working directory to the WaComM++ building directory
3) Create the input, processed, output, and restarts directories
```bash
mkdir input processed output restart results
```
4) Link the source file
```bash
ln -sf ../examples/sources-webinar.json sources.json
```
5) Link the "dry" configuration file
```bash
ln -sf ../examples/webinar-roms-usecase-download.json wacomm.json
```
6) Launch WaComM++ in dry mode
```bash
./wacommplusplus
```
WaComM++ will perform a dry run (the model will not calculate the particles' transport and diffusion)
7) Link the configuration file
```bash
ln -sf ../examples/webinar-native-usecase.json wacomm.json
```
Now WaComM++ is ready to run.

# Running

In our demo, we will run WaComM++ with SLURM manager using 1 Thread and 1 Process:
```bash
sbatch ../examples/slurm_UNITO.sh
```

Watch the output results file:
```bash
tail -f results/wppDemo.out
```

A correct result file will look like this:
```bash
23-12-12 16:05:47,611.037 INFO  WaComMINFO - WaComM - C++ Version
23-12-12 16:05:47,611.482 INFO  WaComMINFO - Parallel: Distributed Memory
23-12-12 16:05:47,611.498 INFO  WaComMINFO - Parallel: Shared Memory
23-12-12 16:05:47,611.506 INFO  WaComMINFO - Acceleration: None
23-12-12 16:05:47,611.521 INFO  WaComMINFO - 0: Using 1/1 processes, each on 112 threads.
23-12-12 16:05:47,611.530 INFO  WaComMINFO - Configuration: wacomm.json
23-12-12 16:05:47,619.909 INFO  WaComMINFO - 0: Input from Ocean Model: processed//ocm3_d03_20190401Z08.nc
23-12-12 16:05:48,182.164 INFO  WaComMINFO - Reading from json:sources.json
23-12-12 16:05:48,401.108 INFO  WaComMINFO - ocean_time:1
23-12-12 16:05:48,401.581 INFO  WaComMINFO - s_rho:30
23-12-12 16:05:48,401.628 INFO  WaComMINFO - eta_rho:1135 xi_rho:1528
23-12-12 16:05:48,453.523 INFO  WaComMINFO - Simulating:1.60491e+09 - 20190401Z08
23-12-12 16:05:48,453.570 INFO  WaComMINFO - Sources:1
23-12-12 16:05:48,531.107 INFO  WaComMINFO - Emitted particles: 250000 Total particles: 250000
23-12-12 16:05:48,578.449 INFO  WaComMINFO - 0: Local particles:250000
23-12-12 16:05:48,902.940 INFO  WaComMINFO - Locally processed 250000 in 0.308113 seconds (811389 particles/second).
23-12-12 16:05:48,919.570 INFO  WaComMINFO - Removed particles: 189774 Particles: 60226
23-12-12 16:05:48,919.586 INFO  WaComMINFO - Evaluate concentration
23-12-12 16:05:48,973.337 INFO  WaComMINFO - EVALUATE sec: 0.0532315
23-12-12 16:05:48,973.401 INFO  WaComMINFO - Processed 250000 in 0.519956 seconds (480810 particles/second).
23-12-12 16:05:48,977.088 INFO  WaComMINFO - Saving restart:restart/WACOMM_rst_20190401Z09
23-12-12 16:05:48,977.307 INFO  WaComMINFO - Preparing NetCDF...
23-12-12 16:05:49,026.408 INFO  WaComMINFO - Saving NetCDF in: restart/WACOMM_rst_20190401Z09.nc
23-12-12 16:05:49,060.070 INFO  WaComMINFO - Saving history:output/wacomm_his_20190401Z08.nc
23-12-12 16:05:49,060.211 INFO  WaComMINFO - Saving in: output/wacomm_his_20190401Z08.nc
23-12-12 16:05:49,067.543 INFO  WaComMINFO - --------------: output/wacomm_his_20190401Z08.nc
23-12-12 16:05:49,187.887 INFO  WaComMINFO - 0: Input from Ocean Model: processed//ocm3_d03_20190401Z09.nc

...
```

1) The introductive section describes used hardware and configuration chosen:
```bash
3-12-12 16:05:47,611.037 INFO  WaComMINFO - WaComM - C++ Version
23-12-12 16:05:47,611.482 INFO  WaComMINFO - Parallel: Distributed Memory
23-12-12 16:05:47,611.498 INFO  WaComMINFO - Parallel: Shared Memory
23-12-12 16:05:47,611.506 INFO  WaComMINFO - Acceleration: None
23-12-12 16:05:47,611.521 INFO  WaComMINFO - 0: Using 1/1 processes, each on 112 threads.

...
```

2) Then, each hour of running can be divided into three sections. The first one shows the configuration used:
```bash
...

23-12-12 16:05:47,619.909 INFO  WaComMINFO - 0: Input from Ocean Model: processed//ocm3_d03_20190401Z08.nc
23-12-12 16:05:48,182.164 INFO  WaComMINFO - Reading from json:sources.json
23-12-12 16:05:48,401.108 INFO  WaComMINFO - ocean_time:1
23-12-12 16:05:48,401.581 INFO  WaComMINFO - s_rho:30
23-12-12 16:05:48,401.628 INFO  WaComMINFO - eta_rho:1135 xi_rho:1528
23-12-12 16:05:48,453.523 INFO  WaComMINFO - Simulating:1.60491e+09 - 20190401Z08
23-12-12 16:05:48,453.570 INFO  WaComMINFO - Sources:1
23-12-12 16:05:48,531.107 INFO  WaComMINFO - Emitted particles: 250000 Total particles: 250000

...
```

3) The second one contains the outer cycle and its results:
```bash
...

23-12-12 16:05:48,578.449 INFO  WaComMINFO - 0: Local particles:250000
23-12-12 16:05:48,902.940 INFO  WaComMINFO - Locally processed 250000 in 0.308113 seconds (811389 particles/second).
23-12-12 16:05:48,919.570 INFO  WaComMINFO - Removed particles: 189774 Particles: 60226
23-12-12 16:05:48,919.586 INFO  WaComMINFO - Evaluate concentration
23-12-12 16:05:48,973.401 INFO  WaComMINFO - Processed 250000 in 0.519956 seconds (480810 particles/second).

...
```

4) The last one saves output NetCDF files and loads the next input file:
```bash
...

23-12-12 16:05:48,977.088 INFO  WaComMINFO - Saving restart:restart/WACOMM_rst_20190401Z09
23-12-12 16:05:48,977.307 INFO  WaComMINFO - Preparing NetCDF...
23-12-12 16:05:49,026.408 INFO  WaComMINFO - Saving NetCDF in: restart/WACOMM_rst_20190401Z09.nc
23-12-12 16:05:49,060.070 INFO  WaComMINFO - Saving history:output/wacomm_his_20190401Z08.nc
23-12-12 16:05:49,060.211 INFO  WaComMINFO - Saving in: output/wacomm_his_20190401Z08.nc
23-12-12 16:05:49,067.543 INFO  WaComMINFO - --------------: output/wacomm_his_20190401Z08.nc
23-12-12 16:05:49,187.887 INFO  WaComMINFO - 0: Input from Ocean Model: processed//ocm3_d03_20190401Z09.nc

...
```