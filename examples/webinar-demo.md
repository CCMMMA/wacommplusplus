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
mkdir input processed output restarts
```
4) Link the sources file
```bash
ln -sf ../examples/sources-webinar.json sources.json
```
5) Link the configuration file
```bash
ln -sf ../examples/webinar-roms-usecase-download.json wacomm.json
```
6) Launch WaComM++
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
Be sure to have a configuration file.
WaComM++ can read Fortran namelists used by classic WaComM implementation (https://github.com/ccmmma/wacomm).
If possible, use native json configuration file.
(WaComM++ automatically searchs for namelist.wacomm or wacomm.json configuration file)


## Shared memory parallelism (OpenMP)
WaComM++ supports shared memory parallelization using OpenMP threads.
Select the number of threads to use exporting the OMP_NUM_THREADS environment variable.
Due to the embarassing parallel algorithm, there is no limitation in the number of used threads.
If the OMP_NUM_THREADS is not specified, OpenMPI assumes the number of threads is equals to the number of available CPU
cores.

```bash
export OMP_NUM_THREADS=n
./wacommplusplus
```

