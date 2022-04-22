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
```bash
ln -sf ../examples/webinar-roms-usecase-download.json wacomm.json
```
(NOT REQUIRED FOR DEMO) Launch WaComM++ in "dry" mode, to download required data.
```bash
./wacommplusplus
```
5) Copy the required files:
```bash
./wacommplusplus
```
WaComM++ will perform a dry run (the model will not actually calculate the particles transport and diffusion).
The demo data files will be downloaded from [meteo@uniparthenope](http://data.meteo.uniparthenope.it:/opendap/opendap/wcm3/d04/),
preprocessed, and saved in processed directory.

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

