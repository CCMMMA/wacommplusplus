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
- USE_OPENACC - Activate the OpenACC acceleration (future feature).
- USE_CUDA - Activate the CUDA acceleration (future feature).
- DEBUG - Add logging printouts (do not use for production or evaluation.)

(WaComM++ uses cmake version 3. In some Linux system cmake is the version 3 by default. In other systems the version 3
must be specified using the cmake3 command. )

- Building on PurpleJeans server requires to set the environement:
```bash
module load gcc-8.3.1 
module load ompi-4.1.0-gcc-8.3.1 
module load cmake/3.19.2 
module load cuda/10.1 
```

- We are going to compile with CUDA, OpenMP and MPI support, but without extra debug logging:
```bash
cmake -DUSE_OMP=ON -DUSE_MPI=ON -DUSE_CUDA=ON -DDEBUG=OFF ..
```

6) Run make and wait

```bash
make
```

# Testing
Download the data.

1) Open a console application (I.e. Terminal)
2) Change the current working directory as the WaComM++ building directory
3) Create the input, processed, output and restarts directories
```bash
mkdir input processed output restarts results
```
4) Link the sources file
```bash
ln -sf ../examples/sources-webinar.json sources.json
```
(NOT REQUIRED FOR DEMO) Link the "dry" configuration file. In this demo the dry configuration will be skipped because needed files have been already downoaded.
WaComM++ will perform a dry run (the model will not actually calculate the particles transport and diffusion).
The demo data files will be downloaded from [meteo@uniparthenope](http://data.meteo.uniparthenope.it:/opendap/opendap/wcm3/d04/),
preprocessed, and saved in processed directory.
```bash
ln -sf ../examples/webinar-roms-usecase-download.json wacomm.json
```
```bash
./wacommplusplus
```

5) Copy the required files:
```bash
cp -r /tmp/processed/ processed/ 
```


7) Link the configuration file
```bash
ln -sf ../examples/webinar-native-usecase.json wacomm.json
```
Now WaComM++ is ready to run.

# Running

In our demo we will run WaComM++ with SLURM manager using 1 Thread, 1 Process and 1 GPU:

```bash
sbatch --gres=gpu:tesla:1 --job-name=wpp1G --output=results/wppDemo.out --error=results/wppDemo.err -c 1 -n 1 ../examples/slurm_webinarGPU.sh
```

Open output results file:
```bash
vi results/wppDemo.out
```

A correct result file will look like:
```bash
22-04-22 11:43:41,405.889 INFO  WaComMINFO - WaComM - C++ Version
22-04-22 11:43:41,406.115 INFO  WaComMINFO - Parallel: Distributed Memory
22-04-22 11:43:41,406.147 INFO  WaComMINFO - Parallel: Shared Memory
22-04-22 11:43:41,406.185 INFO  WaComMINFO - Acceleration: CUDA 1 device(s)
22-04-22 11:43:41,406.477 INFO  WaComMINFO - Device Number: 0
22-04-22 11:43:41,406.510 INFO  WaComMINFO - Device name: Quadro RTX 6000
22-04-22 11:43:41,406.544 INFO  WaComMINFO - Memory Clock Rate (KHz): 7001000
22-04-22 11:43:41,406.568 INFO  WaComMINFO - Memory Bus Width (bits): 384
22-04-22 11:43:41,406.636 INFO  WaComMINFO - Peak Memory Bandwidth (GB/s): 672.096
22-04-22 11:43:41,406.666 INFO  WaComMINFO - 0: Using 1/1 processes, each on 20 threads.
22-04-22 11:43:41,406.702 INFO  WaComMINFO - Configuration: wacomm.json
22-04-22 11:43:41,428.139 INFO  WaComMINFO - 0: Input from Ocean Model: processed//ocm3_d03_20190401Z1000.nc
22-04-22 11:43:43,360.252 INFO  WaComMINFO - Reading from json:sources.json
22-04-22 11:43:44,056.787 INFO  WaComMINFO - ocean_time:1
22-04-22 11:43:44,056.849 INFO  WaComMINFO - s_rho:30
22-04-22 11:43:44,056.891 INFO  WaComMINFO - eta_rho:1135 xi_rho:1528
22-04-22 11:43:44,210.326 INFO  WaComMINFO - Simulating:1.60492e+09 - 20190401Z10
22-04-22 11:43:44,210.490 INFO  WaComMINFO - Sources:1
22-04-22 11:43:44,225.355 INFO  WaComMINFO - Emitted particles: 20000 Total particles: 20000
22-04-22 11:43:44,237.979 INFO  WaComMINFO - 0: Local particles:20000
22-04-22 11:43:45,409.944 INFO  WaComMINFO - Locally processed 20000 in 0.179037 seconds (111709 particles/second).
22-04-22 11:43:45,411.862 INFO  WaComMINFO - Removed particles: 0 Particles: 20000
22-04-22 11:43:45,411.883 INFO  WaComMINFO - Evaluate concentration
22-04-22 11:43:45,418.262 INFO  WaComMINFO - Processed 20000 in 1.20799 seconds (16556.4 particles/second).
22-04-22 11:43:45,418.604 INFO  WaComMINFO - Saving restart:restart/WACOMM_rst_20190401Z11
22-04-22 11:43:51,077.773 INFO  WaComMINFO - Saving history:output/wacomm_his_20190401Z10.nc
22-04-22 11:43:51,077.880 INFO  WaComMINFO - Saving in: output/wacomm_his_20190401Z10.nc
22-04-22 11:43:51,424.383 INFO  WaComMINFO - --------------: output/wacomm_his_20190401Z10.nc
22-04-22 11:43:52,995.869 INFO  WaComMINFO - 0: Input from Ocean Model: processed//ocm3_d03_20190401Z1100.nc

...
```

1)The introductive section describe used hardware and configuration chose:
```bash
22-04-22 11:43:41,405.889 INFO  WaComMINFO - WaComM - C++ Version
22-04-22 11:43:41,406.115 INFO  WaComMINFO - Parallel: Distributed Memory
22-04-22 11:43:41,406.147 INFO  WaComMINFO - Parallel: Shared Memory
22-04-22 11:43:41,406.185 INFO  WaComMINFO - Acceleration: CUDA 1 device(s)
22-04-22 11:43:41,406.477 INFO  WaComMINFO - Device Number: 0
22-04-22 11:43:41,406.510 INFO  WaComMINFO - Device name: Quadro RTX 6000
22-04-22 11:43:41,406.544 INFO  WaComMINFO - Memory Clock Rate (KHz): 7001000
22-04-22 11:43:41,406.568 INFO  WaComMINFO - Memory Bus Width (bits): 384
22-04-22 11:43:41,406.636 INFO  WaComMINFO - Peak Memory Bandwidth (GB/s): 672.096
22-04-22 11:43:41,406.666 INFO  WaComMINFO - 0: Using 1/1 processes, each on 20 threads.

...
```

2) Then each hour of run can be divided into three sections. The first one show the configuration used:
```bash
...

22-04-22 11:43:41,428.139 INFO  WaComMINFO - 0: Input from Ocean Model: processed//ocm3_d03_20190401Z1000.nc
22-04-22 11:43:43,360.252 INFO  WaComMINFO - Reading from json:sources.json
22-04-22 11:43:44,056.787 INFO  WaComMINFO - ocean_time:1
22-04-22 11:43:44,056.849 INFO  WaComMINFO - s_rho:30
22-04-22 11:43:44,056.891 INFO  WaComMINFO - eta_rho:1135 xi_rho:1528
22-04-22 11:43:44,210.326 INFO  WaComMINFO - Simulating:1.60492e+09 - 20190401Z10
22-04-22 11:43:44,210.490 INFO  WaComMINFO - Sources:1
22-04-22 11:43:44,225.355 INFO  WaComMINFO - Emitted particles: 20000 Total particles: 20000

...
```

3) Then each hour of run can be divided into three sections. The second one contains the outer cycle and its results:
```bash
...

22-04-22 11:43:44,237.979 INFO  WaComMINFO - 0: Local particles:20000
22-04-22 11:43:45,409.944 INFO  WaComMINFO - Locally processed 20000 in 0.179037 seconds (111709 particles/second).
22-04-22 11:43:45,411.862 INFO  WaComMINFO - Removed particles: 0 Particles: 20000
22-04-22 11:43:45,411.883 INFO  WaComMINFO - Evaluate concentration
22-04-22 11:43:45,418.262 INFO  WaComMINFO - Processed 20000 in 1.20799 seconds (16556.4 particles/second).

...
```

4) Then each hour of run can be divided into three sections. The last one save output NetCDF files and load the next input file:
```bash
...

22-04-22 11:43:45,418.604 INFO  WaComMINFO - Saving restart:restart/WACOMM_rst_20190401Z11
22-04-22 11:43:51,077.773 INFO  WaComMINFO - Saving history:output/wacomm_his_20190401Z10.nc
22-04-22 11:43:51,077.880 INFO  WaComMINFO - Saving in: output/wacomm_his_20190401Z10.nc
22-04-22 11:43:51,424.383 INFO  WaComMINFO - --------------: output/wacomm_his_20190401Z10.nc
22-04-22 11:43:52,995.869 INFO  WaComMINFO - 0: Input from Ocean Model: processed//ocm3_d03_20190401Z1100.nc

...
```
