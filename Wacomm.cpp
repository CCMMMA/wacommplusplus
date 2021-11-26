//
// Created by Raffaele Montella on 12/5/20.
//

#include "Wacomm.hpp"

#include <chrono>
#include <utility>
#include "OceanModelAdapters/ROMSAdapter.hpp"
#include "JulianDate.hpp"

#ifdef USE_MPI
#define OMPI_SKIP_MPICXX
#include <mpi.h>
#endif

#ifdef USE_OMP
#include <omp.h>
#endif

#ifdef USE_CUDA
#include <cuda.h>
#include <cuda_runtime_api.h>
#include "cuda/kernel.h"

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"CUDA ERROR: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}
#endif

Wacomm::Wacomm(std::shared_ptr<Config> config,
               std::shared_ptr<OceanModelAdapter> oceanModelAdapter,
               std::shared_ptr<Sources> sources,
               std::shared_ptr<Particles> particles) :
               config(std::move(config)),
               oceanModelAdapter(std::move(oceanModelAdapter)),
               sources(std::move(sources)),
               particles(std::move(particles)) {

    logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("WaComM"));
    //cout << "Wacomm::Wacomm oceanModelAdapter->H()(650,550):" << oceanModelAdapter->H()(650,550) << endl;

}

void Wacomm::run()
{
    LOG4CPLUS_DEBUG(logger,"Dry mode:" << config->Dry() );
    int num_gpus=0;
    int world_size=1, world_rank=0;
    int ompMaxThreads=1, ompThreadNum=0;


#ifdef USE_MPI
    // Get the size of the MPI world (the number of available processors)
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the number of the current MPI processor
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
#endif

#ifdef USE_OMP
    // Get the number of threads available for computation
    ompMaxThreads=omp_get_max_threads();
#endif

#ifdef USE_CUDA
    // Get the number of available GPU devices
    cudaError_t err = cudaGetDeviceCount(&num_gpus);

    // Check if a GPU is available
    if (num_gpus <= 0 || err != cudaSuccess) {

        // If a GPU is not present, or a problem occurred, set the number of GPUs as 0;
        num_gpus = 0;
    }
#endif

    // Get the size of the time axis
    size_t ocean_time=oceanModelAdapter->OceanTime().Nx();

    // Get the size of the vertical axis
    size_t s_rho=oceanModelAdapter->SRho().Nx();

    // Get the size of the south-north axis
    size_t eta_rho=oceanModelAdapter->Mask().Nx();

    // Get the size of the west-east axis
    size_t xi_rho=oceanModelAdapter->Mask().Ny();

    // Check if this is the main processor
    if (world_rank==0) {
        // Print some diagnostics
        LOG4CPLUS_INFO(logger, "ocean_time:" + std::to_string(ocean_time));
        LOG4CPLUS_INFO(logger, "s_rho:" + std::to_string(s_rho));
        LOG4CPLUS_INFO(logger, "eta_rho:" + std::to_string(eta_rho) + " xi_rho:" + std::to_string(xi_rho));
    }

    // Define the final concentration matrix
    Array4<float> conc(ocean_time,s_rho,eta_rho,xi_rho,0,-(int)s_rho+1,0,0);

    // Set the concentration matrix to 0
    #pragma omp parallel for collapse(4) default(none) shared(ocean_time, s_rho, eta_rho, xi_rho, conc)
    for (int t=0;t<ocean_time;t++) {
        for (int k=-(int)s_rho+1; k<=0; k++) {
            for (int j=0; j<eta_rho; j++) {
                for (int i=0; i<xi_rho; i++) {
                    conc(t,k,j,i)=0;
                }
            }
        }
    }

    // Define a gregorian calendar
    Calendar cal;

    // Total number of particles (valued only if world_rank==0)
    size_t nParticles = -1;

    // For each element in the time axis
    for (int ocean_time_idx = 0; ocean_time_idx < ocean_time; ocean_time_idx++) {

        // Record start time
        auto start = std::chrono::high_resolution_clock::now();

#ifdef USE_MPI
        // Define a vector of integers hosting the number of particles for each processor
        std::unique_ptr<int[]> send_counts = std::make_unique<int[]>(world_size);

        // Define a vector of integers hosting the displacement o
        std::unique_ptr<int[]> displs = std::make_unique<int[]>(world_size);

        // Define the sending buffer
        std::unique_ptr<struct particle_data[]> sendbuf;

        // Define the receiving buffer
        std::unique_ptr<struct particle_data[]> recvbuf;
#endif
        /*
         * Pointer to the local particles vector:
         * Processor 0: points to the whole particles data structure
         * Processor >0: points to the locala particles data structure
         */
        Particles *pLocalParticles;

        // Check if the processor is the number 0
        if (world_rank == 0) {
            
            // Time in "seconds since 1968-05-23 00:00:00"
            double modJulian = oceanModelAdapter->OceanTime()(ocean_time_idx);

            // Convert time in days based
            modJulian = modJulian / 86400;


            JulianDate::fromModJulian(modJulian, cal);

            LOG4CPLUS_INFO(logger, "Simulating:" << oceanModelAdapter->OceanTime()(ocean_time_idx) << " - "
                                                 << cal.asNCEPdate());
            LOG4CPLUS_INFO(logger, "Sources:" << sources->size());

            // Get the number of particles
            size_t nParticles0 = particles->size();

            // Get the number of sources
            int nSources = sources->size();

            // For each source
            for (int idx = 0; idx < nSources; idx++) {

                // Emit particles
                sources->at(idx).emit(config, particles, JulianDate::toModJulian(cal.asNCEPdate()));
            }

            // Get the total number of particles
            nParticles = particles->size();

            LOG4CPLUS_INFO(logger, "Emitted particles: " << (nParticles-nParticles0) << " Total particles: " << nParticles);

#ifdef USE_MPI
            // Get the number of particles to be processed for each process
            size_t particlesPerProcess = nParticles / world_size;

            // Get the number of spare particles for the process with world_rank==0
            size_t spare = nParticles % world_size;

            // The process with world_rank==0 will calculate extra particles
            send_counts[0] = (int)((particlesPerProcess + spare));

            // The process with world_rank==0 will process the first particles
            displs[0]=0;

            LOG4CPLUS_DEBUG(logger, world_rank << ": send_counts[0]=" << send_counts.get()[0] << " displ[0]=" << displs.get()[0]);

            // Prepare counts and displacement for data distribution
            // For each process...
            for (int i = 1; i < world_size; i++) {

                // Set the number of particles per process
                send_counts.get()[i] = (int)(particlesPerProcess);

                // Set the displacement if terms of particle_data size
                displs.get()[i]=send_counts.get()[0]+particlesPerProcess*(i-1);

                LOG4CPLUS_DEBUG(logger, world_rank << ": send_counts[" << i << "]=" << send_counts.get()[i] << " displ[" << i << "]=" << displs.get()[i]);
            }

            // Prepare a buffer of particle_data
            sendbuf=std::make_unique<struct particle_data[]>(nParticles);

            // For each particle...
            for (int idx=0;idx<nParticles;idx++) {

                // Copy data into the buffer
                sendbuf[idx]=particles->at(idx).data();
            }
#endif
        }

#ifdef USE_MPI
        // Broadcast the number of particles for each processor
        MPI_Bcast(send_counts.get(),world_size,MPI_INT,0,MPI_COMM_WORLD);

        // Get the numeber of particles to process for the current processor
        int particlesToProcess=send_counts.get()[world_rank];

        LOG4CPLUS_DEBUG(logger, world_rank << ":" << " particlesToProcess:" << particlesToProcess);



        // Allocate the receiving buffer
        recvbuf=std::make_unique<particle_data[]>(particlesToProcess);

        // Define a variable that will contain the mpiError
        int mpiError;

        // Define a MPI struct miming particle_data struct

        // Set the number of fields
        constexpr std::size_t num_members = 7;

        // Set the cardinality of each field
        int lengths[num_members] = { 1, 1, 1, 1, 1, 1, 1 };

        // Define an array of MPI int containing the offset of each struct field
        MPI_Aint offsets[num_members] = {
                offsetof(struct particle_data, id),
                offsetof(struct particle_data, k),
                offsetof(struct particle_data, j),
                offsetof(struct particle_data, i),
                offsetof(struct particle_data, health),
                offsetof(struct particle_data, age),
                offsetof(struct particle_data, time)
        };

        // Define an array of MPI data type containing the MPI type of each field
        MPI_Datatype types[num_members] = {
                MPI_UNSIGNED_LONG,
                MPI_DOUBLE,
                MPI_DOUBLE,
                MPI_DOUBLE,
                MPI_DOUBLE,
                MPI_DOUBLE,
                MPI_DOUBLE
        };

        // Define a container for the new MPI data type
        MPI_Datatype mpiParticleData;

        // Create the MPI struct
        MPI_Type_create_struct(num_members, lengths, offsets, types, &mpiParticleData);

        // Add the new datatype
        MPI_Type_commit(&mpiParticleData);

        // Distribute to all processes the send buffer
        mpiError=MPI_Scatterv(sendbuf.get(), send_counts.get(), displs.get(), mpiParticleData,
                     recvbuf.get(), particlesToProcess, mpiParticleData, 0, MPI_COMM_WORLD);

        LOG4CPLUS_DEBUG(logger, world_rank << ": mpiError=" << mpiError);

        // Define a container for the particles that the current processor must process
        Particles localParticles;

        // For each particle to process
        for (int idx=0;idx<particlesToProcess;idx++) {

            // Read the particle from the receiving buffer
            Particle particle(recvbuf[idx]);

            // Add the particle to the local container
            localParticles.push_back(particle);
        }

        // Get the pointer to the local container
        pLocalParticles=&localParticles;
#else
        // Get the pointer to the particles' container
        pLocalParticles = particles.get();
        int particlesToProcess = particles->size();
#endif
        LOG4CPLUS_INFO(logger, world_rank << ": Local particles:" << pLocalParticles->size());

        // Get the number of particles to be processed by each thread
        size_t particlesPerThread = particlesToProcess / ompMaxThreads;

        // Get the number of spare particles for the thread with tidx==0
        size_t sparePerThread = particlesToProcess % ompMaxThreads;

        // Dafine an array with the number of particles to be processed by each thread
        size_t thread_counts[ompMaxThreads];

        // Define an array with the displacement of particles to be processed by each thread
        size_t thread_displs[ompMaxThreads];

        // The first thread get the spare
        thread_counts[0] = (int)((particlesPerThread + sparePerThread));

        // The first tread starts from 0
        thread_displs[0] = 0;

        // For each available thread...
        for (int tidx = 1; tidx < ompMaxThreads; tidx++) {

            // Set the number of particles per thread
            thread_counts[tidx] = (int)(particlesPerThread);

            // Set the displaement
            thread_displs[tidx]=thread_counts[0]+particlesPerThread*(tidx-1);
        }

	    typedef struct WacommVariables{
            double *oceanTimeDevice, *depthIntervalsDevice, *lonRadDevice, *latRadDevice, *maskDevice, *hDevice;
            float *zetaDevice, *uDevice, *vDevice, *wDevice, *aktDevice;
            struct config_data *configDevice;
        } WacommVariables;

	    struct particle_data *particlesHost;
	    WacommVariables *stateVector;

#ifdef USE_CUDA
        // Check if the number of GPUs is more than 0
        if (num_gpus>0) {

            // A pointer to the configuration data
            config_data *pConfigData = config->dataptr();

            // The size of the particle_data struct in bytes
            int itemSize = sizeof(struct particle_data);

            // The size of all particles in byte
            int sizeParticles = pLocalParticles->size() * itemSize;

            // Allocate a container of each GPU
            stateVector = new WacommVariables[num_gpus];

            // Create a host copy of the particles data
            particlesHost = (struct particle_data*)malloc(pLocalParticles->size() * itemSize);

            // Copy each particle
            for(int i=0; i<pLocalParticles->size(); i++){
                particlesHost[i] = pLocalParticles->at(i).data();
            }
            // For each GPU...
            for (int i=0; i<num_gpus; i++){

                //set GPU Device
                cudaSetDevice(i);

                // Get the GPU id
                int gpu_id = -1;
                cudaGetDevice(&gpu_id);

                //device memory allocation
                cudaMalloc((void**) &(stateVector[i].oceanTimeDevice), oceanModelAdapter->OceanTime().Nx() *sizeof(double));
                cudaMalloc((void**) &(stateVector[i].depthIntervalsDevice), oceanModelAdapter->DepthIntervals().Nx() * sizeof(double));
                cudaMalloc((void**) &(stateVector[i].lonRadDevice), oceanModelAdapter->LonRad().Nx() * oceanModelAdapter->LonRad().Ny() * sizeof(double));
                cudaMalloc((void**) &(stateVector[i].latRadDevice), oceanModelAdapter->LatRad().Nx() * oceanModelAdapter->LatRad().Ny() * sizeof(double));
                cudaMalloc((void**) &(stateVector[i].maskDevice), oceanModelAdapter->Mask().Nx() * oceanModelAdapter->Mask().Ny() * sizeof(double));
                cudaMalloc((void**) &(stateVector[i].hDevice), oceanModelAdapter->H().Nx() * oceanModelAdapter->H().Ny() * sizeof(double));
                cudaMalloc((void**) &(stateVector[i].zetaDevice), oceanModelAdapter->Zeta().Nx() * oceanModelAdapter->Zeta().Ny() * oceanModelAdapter->Zeta().Nz() * sizeof(float));
                cudaMalloc((void**) &(stateVector[i].uDevice), oceanModelAdapter->U().Nx() * oceanModelAdapter->U().Ny() * oceanModelAdapter->U().Nz() * oceanModelAdapter->U().N4() * sizeof(float));
                cudaMalloc((void**) &(stateVector[i].vDevice), oceanModelAdapter->V().Nx() * oceanModelAdapter->V().Ny() * oceanModelAdapter->V().Nz() * oceanModelAdapter->V().N4() * sizeof(float));
                cudaMalloc((void**) &(stateVector[i].wDevice), oceanModelAdapter->W().Nx() * oceanModelAdapter->W().Ny() * oceanModelAdapter->W().Nz() * oceanModelAdapter->W().N4() * sizeof(float));
                cudaMalloc((void**) &(stateVector[i].aktDevice), oceanModelAdapter->AKT().Nx() * oceanModelAdapter->AKT().Ny() * oceanModelAdapter->AKT().Nz() * oceanModelAdapter->AKT().N4() * sizeof(float));
                cudaMalloc((void**) &(stateVector[i].configDevice), sizeof(struct config_data));

                //copy data from host to device
                cudaMemcpy(stateVector[i].oceanTimeDevice, oceanModelAdapter->OceanTime(), oceanModelAdapter->OceanTime().Nx() * sizeof(double), cudaMemcpyHostToDevice);
                cudaMemcpy(stateVector[i].depthIntervalsDevice, oceanModelAdapter->DepthIntervals(), oceanModelAdapter->DepthIntervals().Nx()*sizeof(double), cudaMemcpyHostToDevice);
                cudaMemcpy(stateVector[i].lonRadDevice, oceanModelAdapter->LonRad(), oceanModelAdapter->LonRad().Nx() * oceanModelAdapter->LonRad().Ny() * sizeof(double), cudaMemcpyHostToDevice);
                cudaMemcpy(stateVector[i].latRadDevice, oceanModelAdapter->LatRad(), oceanModelAdapter->LatRad().Nx() * oceanModelAdapter->LatRad().Ny() * sizeof(double), cudaMemcpyHostToDevice);
                cudaMemcpy(stateVector[i].maskDevice, oceanModelAdapter->Mask(), oceanModelAdapter->Mask().Nx() * oceanModelAdapter->Mask().Ny() * sizeof(double), cudaMemcpyHostToDevice);
                cudaMemcpy(stateVector[i].hDevice, oceanModelAdapter->H(), oceanModelAdapter->H().Nx() * oceanModelAdapter->H().Ny() * sizeof(double), cudaMemcpyHostToDevice);
                cudaMemcpy(stateVector[i].zetaDevice, oceanModelAdapter->Zeta(), oceanModelAdapter->Zeta().Nx() * oceanModelAdapter->Zeta().Ny() * oceanModelAdapter->Zeta().Nz() * sizeof(float), cudaMemcpyHostToDevice);
                cudaMemcpy(stateVector[i].uDevice, oceanModelAdapter->U(),
                           oceanModelAdapter->U().Nx() * oceanModelAdapter->U().Ny() * oceanModelAdapter->U().Nz() *
                           oceanModelAdapter->U().N4() *
                           sizeof(float), cudaMemcpyHostToDevice);
                cudaMemcpy(stateVector[i].vDevice, oceanModelAdapter->V(),
                           oceanModelAdapter->V().Nx() * oceanModelAdapter->V().Ny() * oceanModelAdapter->V().Nz() *
                           oceanModelAdapter->V().N4() *
                           sizeof(float), cudaMemcpyHostToDevice);
                cudaMemcpy(stateVector[i].wDevice, oceanModelAdapter->W(),
                           oceanModelAdapter->W().Nx() * oceanModelAdapter->W().Ny() * oceanModelAdapter->W().Nz() *
                           oceanModelAdapter->W().N4() *
                           sizeof(float), cudaMemcpyHostToDevice);
                cudaMemcpy(stateVector[i].aktDevice, oceanModelAdapter->AKT(),
                           oceanModelAdapter->AKT().Nx() * oceanModelAdapter->AKT().Ny() * oceanModelAdapter->AKT().Nz() *
                           oceanModelAdapter->AKT().N4() *
                           sizeof(float), cudaMemcpyHostToDevice);
                cudaMemcpy(stateVector[i].configDevice, pConfigData, sizeof(struct config_data), cudaMemcpyHostToDevice);
            }
        }
#endif
        // Record start time
        auto startLocal = std::chrono::high_resolution_clock::now();
        // Begin the shared memory parallel section
        #pragma omp parallel default(none) private(ompThreadNum) shared(thread_counts, thread_displs, config, pLocalParticles, particlesPerThread, ocean_time_idx, oceanModelAdapter, num_gpus, particlesHost, stateVector)
        {

#ifdef USE_OMP
            // Get the number of the current thread
            ompThreadNum = omp_get_thread_num();
#endif

            // Get the index (array pLocalParticles) of the first particle the thread must process
            size_t first = thread_displs[ompThreadNum];

            // Get the index (array pLocalParticles) of the last particle the thread must process
            size_t last = first + thread_counts[ompThreadNum];

            // Check if no GPU is available
            if (num_gpus<=0) {

		        // For each particle the thread must process...
            	for (size_t idx = first; idx < last; idx++) {

                	// Move the particle
                	pLocalParticles->at(idx).move(config->dataptr(),
                      	ocean_time_idx,
                      	oceanModelAdapter->OceanTime(),
                      	oceanModelAdapter->Mask(),
                      	oceanModelAdapter->LonRad(),
                      	oceanModelAdapter->LatRad(),
                      	oceanModelAdapter->SW(),
                      	oceanModelAdapter->DepthIntervals(),
                      	oceanModelAdapter->H(),
                      	oceanModelAdapter->Zeta(),
                      	oceanModelAdapter->U(),
                      	oceanModelAdapter->V(),
                      	oceanModelAdapter->W(),
                      	oceanModelAdapter->AKT()
                  );
		        }
            }
#ifdef USE_CUDA
	        else {
                if (num_gpus>0) {
                    int particlesToProcessGPU = particlesPerThread;

                    if (particlesToProcessGPU > 0){
                        // Get the number of particles to be processed by each GPU
                        size_t particlesPerGPU = particlesToProcessGPU / num_gpus;

                        // Get the number of spare particles for the GPU with idx==0
                        size_t sparePerGPU = particlesToProcessGPU % num_gpus;

                        // Dafine an array with the number of particles to be processed by each GPU
                        size_t GPU_counts[num_gpus];

                        // Define an array with the displacement of particles to be processed by each GPU
                        size_t GPU_displs[num_gpus];

                        // The first GPU get the spare
                        GPU_counts[0] = (int)((particlesPerGPU + sparePerGPU));

                        // The first GPU starts from 0
                        GPU_displs[0] = 0;

                        // For each available GPU...
                        for (int gidx = 1; gidx < num_gpus; gidx++) {
                            // Set the number of particles per GPU
                            GPU_counts[gidx] = (int)(particlesPerGPU);

                            // Set the displaement
                            GPU_displs[gidx] = GPU_counts[0] + particlesPerGPU*(gidx-1);
                        }

                        typedef struct ThreadSectionDevice{
                            struct particle_data *sectionParticlesDevice;
                        } ThreadSectionDevice;

                        ThreadSectionDevice *threadSectionDevice = new ThreadSectionDevice[num_gpus];

                        for (int idx=0; idx < num_gpus; idx++){
                            size_t GPU_first = GPU_displs[idx] + first;
                            size_t GPU_last = GPU_first + GPU_counts[idx];

                            //printf("GPU_first: %d, GPU_last: %d\n", GPU_first, GPU_last);

                            struct particle_data *particlesThread = &particlesHost[GPU_first];

                            int gpu_id = -1;
                            cudaSetDevice(idx);
                            cudaGetDevice(&gpu_id);

                            cudaMalloc((void**) &(threadSectionDevice[gpu_id].sectionParticlesDevice), GPU_counts[idx] * sizeof(struct particle_data));
                            cudaMemcpy(threadSectionDevice[gpu_id].sectionParticlesDevice, particlesThread, GPU_counts[idx] * sizeof(struct particle_data), cudaMemcpyHostToDevice);

                            gpuErrchk(cudaMoveParticle(stateVector[gpu_id].configDevice,
                                threadSectionDevice[gpu_id].sectionParticlesDevice,
                                ocean_time_idx,
                                oceanModelAdapter->OceanTime().Nx(),
                                oceanModelAdapter->SW().Nx(),
                                oceanModelAdapter->SRho().Nx(),
                                oceanModelAdapter->Mask().Nx(),
                                oceanModelAdapter->Mask().Ny(),
                                stateVector[gpu_id].oceanTimeDevice,
                                stateVector[gpu_id].maskDevice,
                                stateVector[gpu_id].lonRadDevice,
                                stateVector[gpu_id].latRadDevice,
                                stateVector[gpu_id].depthIntervalsDevice,
                                stateVector[gpu_id].hDevice,
                                stateVector[gpu_id].zetaDevice,
                                stateVector[gpu_id].uDevice,
                                stateVector[gpu_id].vDevice,
                                stateVector[gpu_id].wDevice,
                                stateVector[gpu_id].aktDevice,
                                GPU_counts[idx], idx, gpu_id));

                            //copy from device to host
                            cudaMemcpy(particlesThread, threadSectionDevice[gpu_id].sectionParticlesDevice, GPU_counts[idx] * sizeof(struct particle_data), cudaMemcpyDeviceToHost);

                            // Copy all thread particles to the local processor particles
                            for(int i=GPU_first; i < GPU_last; i++){
                                //printf("health: %f\n", particlesThread[k].health);
                                pLocalParticles->at(i).data(particlesThread[i - GPU_first]);
                            }

                            cudaFree(threadSectionDevice[gpu_id].sectionParticlesDevice);
                        }
                    }
                }
            }
#endif
        } 
            
        // Record end time
        auto finishLocal = std::chrono::high_resolution_clock::now();

#ifdef USE_MPI
        // For each particle processed by the current process...
        for (int idx=0;idx<localParticles.size();idx++) {
            // Copy the particle data in the receivig buffer
            recvbuf[idx]=localParticles.at(idx).data();
        }

        // Send the receiving buffer to the process with world_rank==0
        MPI_Gatherv(recvbuf.get(), send_counts.get()[world_rank],mpiParticleData,
                    sendbuf.get(),send_counts.get(),displs.get(),mpiParticleData,0,MPI_COMM_WORLD);

        // Remove the MPI Data type
        MPI_Type_free(&mpiParticleData);
#endif

#ifdef USE_CUDA
        if (num_gpus>0) {
            free(particlesHost);

            for(int i=0; i<num_gpus; i++){
                cudaSetDevice(i);

                cudaFree(stateVector[i].depthIntervalsDevice);
                cudaFree(stateVector[i].oceanTimeDevice);
                cudaFree(stateVector[i].lonRadDevice);
                cudaFree(stateVector[i].latRadDevice);
                cudaFree(stateVector[i].maskDevice);
                cudaFree(stateVector[i].hDevice);
                cudaFree(stateVector[i].zetaDevice);
                cudaFree(stateVector[i].uDevice);
                cudaFree(stateVector[i].vDevice);
                cudaFree(stateVector[i].wDevice);
                cudaFree(stateVector[i].aktDevice);
                cudaFree(stateVector[i].configDevice);
            }
        }
#endif
        // Check if the current process is the world_rank==0
        if (world_rank==0) {
            // Calculate the local wall clock
            std::chrono::duration<double> elapsedLocal = finishLocal - startLocal;

            // Calculate the local particles per second
            double nParticlesPerSecondLocal = pLocalParticles->size() / elapsedLocal.count();
            LOG4CPLUS_INFO(logger, "Locally processed " << pLocalParticles->size() << " in " << elapsedLocal.count()
                                                        << " seconds (" << nParticlesPerSecondLocal
                                                        << " particles/second).");

#ifdef USE_MPI
            // For each particle...
            for (int idx = 0; idx < particles->size(); idx++) {
                // Copy the particle status from the sending buffer
                particles->at(idx).data(sendbuf[idx]);
            }
#endif
            // Get the old number of particles
            size_t nParticles0 = particles->size();

            // Remove dead particles
            particles->erase(std::remove_if(particles->begin(), particles->end(),
                                            [](const Particle &particle) { return !particle.isAlive(); }),
                             particles->end());

            size_t nParticles = particles->size();
            LOG4CPLUS_INFO(logger, "Removed particles: " << (nParticles0 - nParticles) << " Particles: " << nParticles);

            LOG4CPLUS_INFO(logger, "Evaluate concentration");

            // Evaluate the concentration of particles per grid cell
#pragma omp parallel for default(none) shared(nParticles, ocean_time_idx, s_rho, eta_rho, xi_rho, conc)
            // For each particle...
            for (int idx = 0; idx < nParticles; idx++) {
                // Get the reference to the particle
                const Particle &particle = particles->at(idx);

                // Even if it is redundant, check if the particle is alive
                if (particle.isAlive()) {

                    // Get the integer indeces k, j, i
                    int k = (int) round(particle.K());
                    int j = (int) round(particle.J());
                    int i = (int) round(particle.I());

                    // Check if the indices are consistent
                    if (j >= 0 && j < eta_rho && i > 0 && i < xi_rho && k >= (-(int) s_rho + 1) && k <= 0) {

                        // Increment the count of the particles in the grid cell
                        conc(ocean_time_idx, k, j, i) = conc(ocean_time_idx, k, j, i) + 1;
                    }
                }
            }

            if (config->MaskOutput()) {
                // Mask all grid cells belonging to the land
#pragma omp parallel for collapse(3) default(none) shared(ocean_time_idx, s_rho, eta_rho, xi_rho, conc)
                // For each level...
                for (int k = -(int) s_rho + 1; k <= 0; k++) {

                    // For each row...
                    for (int j = 0; j < eta_rho; j++) {

                        // For each column...
                        for (int i = 0; i < xi_rho; i++) {

                            // Check if the cell is land
                            if (oceanModelAdapter->Mask()(j, i) != 1) {

                                // Mask the cell
                                conc(ocean_time_idx, k, j, i) = 1e37;
                            }
                        }
                    }
                }
            }
        }

        // Check if this is the process with world_rank==0
        if (world_rank==0) {

            // Record the process with world_rank==0 end time
            auto finish = std::chrono::high_resolution_clock::now();

            // Calculate the wall clock
            std::chrono::duration<double> elapsed = finish - start;

            // Calculate the number of particles per second
            double nParticlesPerSecond = nParticles / elapsed.count();

            LOG4CPLUS_INFO(logger, "Processed " << nParticles << " in " << elapsed.count() << " seconds ("<< nParticlesPerSecond <<" particles/second).");
        }
    }

    // Check if this is the process with world_rank==0
    if (world_rank==0) {
        Calendar calFinal;

        // Calculate the oceanTime at the end of calculations
        double finalOceanTime=oceanModelAdapter->OceanTime()[ocean_time-1]+config->Deltat();

        // Convert the ocean time from modified julian to gregorian calendar
        JulianDate::fromModJulian(finalOceanTime/86400, calFinal);

        // Checl if the history must be saved
        if (!config->SaveHistory().empty()) {

            // Create the history filename
            string historyFilename = config->HistoryRoot() + calFinal.asNCEPdate() ;

            LOG4CPLUS_INFO(logger, "Saving restart:" << historyFilename);

            // Check if the history has to be saved as text (Fortran WaComM compatibility)
            if (config->SaveHistory() == "text") {

                // Save the history as text
                particles->saveAsTxt(historyFilename+ ".txt");

                // Check if the history has to be saved as geojson - to be used only for a small amount of particles
            } else if (config->SaveHistory() == "json") {

                // Save the history as geojson
                particles->saveAsJson(historyFilename+ ".json", finalOceanTime, oceanModelAdapter);

                // Check if the history has to be saved as NetCDF file (new style, suggested)
            } else if (config->SaveHistory() == "nc") {

                // Save the history as NetCDF
                particles->saveAsNetCDF(historyFilename + ".nc", finalOceanTime, oceanModelAdapter);
            }
        }

        // Create the output filename
        string ncOutputFilename=config->NcOutputRoot()+cal.asNCEPdate()+".nc";

        LOG4CPLUS_INFO(logger, "Saving history:" << ncOutputFilename);

        // Save the history
        save(ncOutputFilename,conc );
        //save(ncOutputFilename);
    }
};

Wacomm::~Wacomm() = default;

void Wacomm::save(const string &fileName, Array4<float> &conc) {

    size_t ocean_time=oceanModelAdapter->OceanTime().Nx();
    size_t s_rho=oceanModelAdapter->SRho().Nx();
    size_t eta_rho=oceanModelAdapter->Mask().Nx();
    size_t xi_rho=oceanModelAdapter->Mask().Ny();



    LOG4CPLUS_INFO(logger,"Saving in: " << fileName);

    // Open the file for read access
    netCDF::NcFile dataFile(fileName, netCDF::NcFile::replace,NcFile::nc4);
    LOG4CPLUS_INFO(logger,"--------------: " << fileName);

    NcDim oceanTimeDim = dataFile.addDim("ocean_time", ocean_time);
    NcVar oceanTimeVar = dataFile.addVar("ocean_time", ncDouble, oceanTimeDim);
    oceanTimeVar.putAtt("long_name","time since initialization");
    oceanTimeVar.putAtt("units","seconds since 1968-05-23 00:00:00 GMT");
    oceanTimeVar.putAtt("calendar","gregorian");
    oceanTimeVar.putAtt("field","time, scalar, series");
    oceanTimeVar.putAtt("_CoordinateAxisType","Time");
    oceanTimeVar.putVar(oceanModelAdapter->OceanTime()());

    NcDim sRhoDim = dataFile.addDim("s_rho", s_rho);
    NcVar sRhoVar = dataFile.addVar("s_rho", ncDouble, sRhoDim);
    sRhoVar.putAtt("long_name","S-coordinate at RHO-points");
    sRhoVar.putAtt("positive","up");
    sRhoVar.putAtt("standard_name","ocean_sigma_coordinates");
    sRhoVar.putAtt("field","s_rho, scalar");
    sRhoVar.putAtt("_CoordinateTransformType","Vertical");
    sRhoVar.putAtt("_CoordinateAxes","s_rho");
    sRhoVar.putAtt("_CoordinateAxisType","GeoZ");
    sRhoVar.putAtt("_CoordinateZisPositive","up");
    sRhoVar.putVar(oceanModelAdapter->SRho()());

    NcDim etaRhoDim = dataFile.addDim("eta_rho", eta_rho);
    NcDim xiRhoDim = dataFile.addDim("eta_xi", xi_rho);
    vector<NcDim> etaRhoXiRhoDims;
    etaRhoXiRhoDims.push_back(etaRhoDim);
    etaRhoXiRhoDims.push_back(xiRhoDim);
    NcVar latRhoVar = dataFile.addVar("lat_rho", ncDouble, etaRhoXiRhoDims);
    latRhoVar.putAtt("long_name","latitude of rho-points");
    latRhoVar.putAtt("unit","degree_north");
    latRhoVar.putAtt("standard_name","latitude");
    latRhoVar.putAtt("field","lat_rho, scalar");
    latRhoVar.putAtt("_coordinateaxistype","lat");
    latRhoVar.putVar(oceanModelAdapter->Lat()());

    NcVar lonRhoVar = dataFile.addVar("lon_rho", ncDouble, etaRhoXiRhoDims);
    lonRhoVar.putAtt("long_name","longitude of rho-points");
    lonRhoVar.putAtt("unit","degree_east");
    lonRhoVar.putAtt("standard_name","longitude");
    lonRhoVar.putAtt("field","lon_rho, scalar");
    lonRhoVar.putAtt("_coordinateaxistype","lon");
    lonRhoVar.putVar(oceanModelAdapter->Lon()());

    NcVar maskRhoVar = dataFile.addVar("mask_rho", ncDouble, etaRhoXiRhoDims);
    maskRhoVar.putAtt("long_name","mask on RHO-points");
    maskRhoVar.putAtt("coordinates","lon_rho lat_rho");
    maskRhoVar.putAtt("units","1");
    maskRhoVar.putVar(oceanModelAdapter->Mask()());

    NcVar hVar = dataFile.addVar("h", ncDouble, etaRhoXiRhoDims);
    hVar.putAtt("long_name","bathymetry at RHO-point");
    hVar.putAtt("units","meter");
    hVar.putAtt("coordinates","lon_rho lat_rho");
    hVar.putAtt("field","bath, scalar");
    hVar.putVar(oceanModelAdapter->H()());

    vector<NcDim> oceanTimeEtaRhoXiRhoDims;
    oceanTimeEtaRhoXiRhoDims.push_back(oceanTimeDim);
    oceanTimeEtaRhoXiRhoDims.push_back(etaRhoDim);
    oceanTimeEtaRhoXiRhoDims.push_back(xiRhoDim);

    NcVar hZeta = dataFile.addVar("zeta", ncFloat, oceanTimeEtaRhoXiRhoDims);
    hZeta.putAtt("long_name","free-surface");
    hZeta.putAtt("units","meter");
    hZeta.putAtt("time","ocean_time");
    hZeta.putAtt("coordinates","lon_rho lat_rho ocean_time");
    hZeta.putAtt("field","free-surface, scalar, series");
    hZeta.putAtt("_FillValue",ncFloat, 9.99999993e+36);
    hZeta.putVar(oceanModelAdapter->Zeta()());


    vector<NcDim> oceanTimeSRhoEtaRhoXiRhoDims;
    oceanTimeSRhoEtaRhoXiRhoDims.push_back(oceanTimeDim);
    oceanTimeSRhoEtaRhoXiRhoDims.push_back(sRhoDim);
    oceanTimeSRhoEtaRhoXiRhoDims.push_back(etaRhoDim);
    oceanTimeSRhoEtaRhoXiRhoDims.push_back(xiRhoDim);

    NcVar concVar = dataFile.addVar("conc", ncFloat, oceanTimeSRhoEtaRhoXiRhoDims);
    concVar.putAtt("long_name","concentration_of_suspended_matter_in_sea_water");
    concVar.putAtt("units","1");
    concVar.putAtt("coordinates","lon_rho lat_rho s_rho ocean_time");
    concVar.putAtt("field","");
    concVar.putAtt("time","ocean_time");
    concVar.putAtt("_FillValue",ncFloat, 9.99999993e+36);
    concVar.putVar(conc());

}
/*
void Wacomm::save(const string &fileName)
{
    size_t ocean_time=oceanModelAdapter->OceanTime().Nx();
    Array4<float> conc(ocean_time,
                       oceanModelAdapter->Depth().Nx(),
                       oceanModelAdapter->Latitude().Nx(),
                       oceanModelAdapter->Longitude().Nx());

    for (int k=0; k<oceanModelAdapter->Depth().Nx(); k++) {
        for (int j=0; k<oceanModelAdapter->Latitude().Nx(); j++) {
            for (int i=0; k<oceanModelAdapter->Longitude().Nx(); i++) {
                conc(0,k,j,i)=0;
            }
        }
    }

    LOG4CPLUS_INFO(logger,"Preparing NetCDF...");

    double maxDepth=oceanModelAdapter->Depth()[0];
    double minDepth=oceanModelAdapter->Depth()[oceanModelAdapter->Depth().Nx()-1];
    double minLat=oceanModelAdapter->Latitude()[0];
    double maxLat=oceanModelAdapter->Latitude()[oceanModelAdapter->Latitude().Nx()-1];
    double minLon=oceanModelAdapter->Longitude()[0];
    double maxLon=oceanModelAdapter->Longitude()[oceanModelAdapter->Longitude().Nx()-1];

    // Get the number of particles
    size_t nParticles = particles->size();

    // For each particle...
    for (int idx = 0; idx < nParticles; idx++) {
        // Get the reference to the particle
        const Particle &particle = particles->at(idx);


        if (particle.J() < 0 || particle.I() < 0 || particle.K() > 0) {
            continue;
        }


        double depth, latitude, longitude;
        oceanModelAdapter->kji2deplatlon(particle.K(), particle.J(), particle.I(), depth, latitude, longitude);

        if (depth == 1e37 || latitude == 1e37 || longitude == 1e37) {
            continue;
        }

        double concK=oceanModelAdapter->Depth().Nx()*depth/(maxDepth-minDepth);
        double concJ=oceanModelAdapter->Latitude().Nx()*latitude/(maxLat-minLat);
        double concI=oceanModelAdapter->Longitude().Nx()*longitude/(maxLon-minLon);

        conc(ocean_time,concK,concJ,concI )++;
    }


    LOG4CPLUS_INFO(logger,"Saving NetCDF in: " << fileName);

    // Open the file for read access
    netCDF::NcFile dataFile(fileName, NcFile::write,NcFile::nc4);

    NcDim oceanTimeDim = dataFile.addDim("ocean_time", ocean_time);
    NcVar oceanTimeVar = dataFile.addVar("ocean_time", ncDouble, oceanTimeDim);
    oceanTimeVar.putAtt("long_name","time since initialization");
    oceanTimeVar.putAtt("units","seconds since 1968-05-23 00:00:00 GMT");
    oceanTimeVar.putAtt("calendar","gregorian");
    oceanTimeVar.putAtt("field","time, scalar, series");
    oceanTimeVar.putAtt("_CoordinateAxisType","Time");
    oceanTimeVar.putVar(oceanModelAdapter->OceanTime()());

    NcDim latitudeDim = dataFile.addDim("latitude", oceanModelAdapter->Latitude().Nx());
    NcVar laitudeVar = dataFile.addVar("latitude", ncDouble, latitudeDim);
    laitudeVar.putAtt("long_name","latitude o");
    laitudeVar.putAtt("unit","degree_north");
    laitudeVar.putAtt("standard_name","latitude");
    laitudeVar.putAtt("field","latitude, scalar");
    laitudeVar.putAtt("_coordinateaxistype","lat");
    laitudeVar.putVar(oceanModelAdapter->Latitude()());

    NcDim longitudeDim = dataFile.addDim("longitude", oceanModelAdapter->Longitude().Nx());
    NcVar longitudeVar = dataFile.addVar("longitude", ncDouble, longitudeDim);
    longitudeVar.putAtt("long_name","longitude o");
    longitudeVar.putAtt("unit","degree_north");
    longitudeVar.putAtt("standard_name","longitude");
    longitudeVar.putAtt("field","latitude, scalar");
    longitudeVar.putAtt("_coordinateaxistype","lon");
    longitudeVar.putVar(oceanModelAdapter->Longitude()());

    vector<NcDim> oceanTimeDepthLatitudeLongitudeDims;
    oceanTimeDepthLatitudeLongitudeDims.push_back(oceanTimeDim);
    oceanTimeDepthLatitudeLongitudeDims.push_back(depthDim);
    oceanTimeDepthLatitudeLongitudeDims.push_back(latitudeDim);
    oceanTimeDepthLatitudeLongitudeDims.push_back(longitudeDim);

    NcVar concVar = dataFile.addVar("conc", ncFloat, oceanTimeDepthLatitudeLongitudeDims);
    concVar.putAtt("long_name","concentration_of_suspended_matter_in_sea_water");
    concVar.putAtt("units","1");
    concVar.putAtt("coordinates","longitude latitude depth ocean_time");
    concVar.putAtt("field","");
    concVar.putAtt("time","ocean_time");
    concVar.putAtt("_FillValue",ncFloat, 9.99999993e+36);
    concVar.putVar(conc());

}
*/
