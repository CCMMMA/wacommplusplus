//
// Created by Raffaele Montella on 12/5/20.
//

#include "Wacomm.hpp"

#include <chrono>
#include <utility>
#include "OceanModelAdapters/ROMSAdapter.hpp"
#include "JulianDate.hpp"

#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef USE_OMP
#include <omp.h>
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
    int world_size=1, world_rank=0;
    int ompMaxThreads=1, ompThreadNum=0;

#ifdef USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
#endif

#ifdef USE_OMP
    ompMaxThreads=omp_get_max_threads();
#endif


    size_t ocean_time=oceanModelAdapter->OceanTime().Nx();
    size_t s_rho=oceanModelAdapter->SRho().Nx();
    size_t eta_rho=oceanModelAdapter->Mask().Nx();
    size_t xi_rho=oceanModelAdapter->Mask().Ny();


    if (world_rank==0) {
        LOG4CPLUS_INFO(logger, "ocean_time:" + std::to_string(ocean_time));
        LOG4CPLUS_INFO(logger, "s_rho:" + std::to_string(s_rho));
        LOG4CPLUS_INFO(logger, "eta_rho:" + std::to_string(eta_rho) + " xi_rho:" + std::to_string(xi_rho));
    }

    Array4<float> conc(ocean_time,s_rho,eta_rho,xi_rho,0,-(int)s_rho+1,0,0);

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

    if (!config->Dry()) {
        int itemSize = sizeof(struct particle_data);
        Calendar cal;



        for (int ocean_time_idx = 0; ocean_time_idx < ocean_time; ocean_time_idx++) {

            // Record start time
            auto start = std::chrono::high_resolution_clock::now();


#ifdef USE_MPI
            std::unique_ptr<int[]> send_counts = std::make_unique<int[]>(world_size);
            std::unique_ptr<int[]> displs = std::make_unique<int[]>(world_size);

            std::unique_ptr<struct particle_data[]> sendbuf;
            std::unique_ptr<struct particle_data[]> recvbuf;
#endif
            Particles *pLocalParticles;

            if (world_rank == 0) {


                // Time in "seconds since 1968-05-23 00:00:00"
                double modJulian = oceanModelAdapter->OceanTime()(ocean_time_idx);

                // Convert time in days based
                modJulian = modJulian / 86400;

                JulianDate::fromModJulian(modJulian, cal);

                LOG4CPLUS_INFO(logger, "Simulating:" << oceanModelAdapter->OceanTime()(ocean_time_idx) << " - "
                                                     << cal.asNCEPdate());
                LOG4CPLUS_INFO(logger, "Sources:" << sources->size());

                // Get the number of sources
                int nSources = sources->size();

                for (int idx = 0; idx < nSources; idx++) {
                    Source &source = sources->at(idx);
                    source.emit(config, particles, oceanModelAdapter->OceanTime()(ocean_time_idx));
                }

                LOG4CPLUS_INFO(logger, "Total particles:" << particles->size());

#ifdef USE_MPI
                size_t elementsPerProcess = particles->size() / world_size;
                size_t spare = particles->size() % world_size;

                send_counts[0] = (int)((elementsPerProcess + spare));
                displs[0]=0;
                LOG4CPLUS_DEBUG(logger, world_rank << ": send_counts[0]=" << send_counts.get()[0] << " displ[0]=" << displs.get()[0]);

                for (int i = 1; i < world_size; i++) {
                    send_counts.get()[i] = (int)(elementsPerProcess);//*itemSize);
                    displs.get()[i]=send_counts.get()[0]*itemSize+elementsPerProcess*(i-1)*itemSize;

                    LOG4CPLUS_DEBUG(logger, world_rank << ": send_counts[" << i << "]=" << send_counts.get()[i] << " displ[" << i << "]=" << displs.get()[i]);
                }

                sendbuf=std::make_unique<struct particle_data[]>(particles->size()*itemSize);
                for (int idx=0;idx<particles->size();idx++) {
                    sendbuf[idx]=particles->at(idx).data();
                }
#endif
            }
            int nParticles = particles->size();
#ifdef USE_MPI
            MPI_Bcast(send_counts.get(),world_size,MPI_INT,0,MPI_COMM_WORLD);
            int elementToProcess=send_counts.get()[world_rank];
            LOG4CPLUS_DEBUG(logger, world_rank << ":" << " elementToProcess:" << elementToProcess);

            recvbuf=std::make_unique<particle_data[]>(elementToProcess*itemSize);

            int mpiError;
            constexpr std::size_t num_members = 7;
            int lengths[num_members] = { 1, 1, 1, 1, 1, 1, 1 };
            MPI_Aint offsets[num_members] = {
                    offsetof(struct particle_data, id),
                    offsetof(struct particle_data, k),
                    offsetof(struct particle_data, j),
                    offsetof(struct particle_data, i),
                    offsetof(struct particle_data, health),
                    offsetof(struct particle_data, age),
                    offsetof(struct particle_data, time)
            };
            MPI_Datatype types[num_members] = { MPI_UNSIGNED_LONG, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE };

            MPI_Datatype mpiParticleData;
            MPI_Type_create_struct(num_members, lengths, offsets, types, &mpiParticleData);
            MPI_Type_commit(&mpiParticleData);
            mpiError=MPI_Scatterv(sendbuf.get(), send_counts.get(), displs.get(), mpiParticleData,
                         recvbuf.get(), elementToProcess, mpiParticleData, 0, MPI_COMM_WORLD);
            LOG4CPLUS_DEBUG(logger, world_rank << ": mpiError=" << mpiError);

            Particles localParticles;
            for (int idx=0;idx<elementToProcess;idx++) {
                Particle particle(recvbuf[idx]);
                localParticles.push_back(particle);
            }

            pLocalParticles=&localParticles;
#else
            pLocalParticles = particles.get();
#endif
            LOG4CPLUS_DEBUG(logger, world_rank << ": Local particles:" << pLocalParticles->size());

            config_data *pConfigData = config->dataptr();
            cout << "Wacomm::run oceanModelAdapter->OceanTime()("<< ocean_time_idx <<"): " << oceanModelAdapter->OceanTime()(ocean_time_idx) << endl;


            Array1<double> oceanTime(oceanModelAdapter->OceanTime().Nx(),
                                     oceanModelAdapter->OceanTime().Ox());

            Array2<double> mask(oceanModelAdapter->Mask().Nx(), oceanModelAdapter->Mask().Ny(),
                                oceanModelAdapter->Mask().Ox(), oceanModelAdapter->Mask().Oy());

            Array2<double> lonRad(oceanModelAdapter->LonRad().Nx(), oceanModelAdapter->LonRad().Ny(),
                                  oceanModelAdapter->LonRad().Ox(), oceanModelAdapter->LonRad().Oy());

            Array2<double> latRad(oceanModelAdapter->LatRad().Nx(), oceanModelAdapter->LatRad().Ny(),
                                  oceanModelAdapter->LatRad().Ox(), oceanModelAdapter->LatRad().Oy());

            Array1<double> depth(oceanModelAdapter->Depth().Nx(),
                                 oceanModelAdapter->Depth().Ox());

            Array2<double> h(oceanModelAdapter->H().Nx(), oceanModelAdapter->H().Ny(),
                             oceanModelAdapter->H().Ox(), oceanModelAdapter->H().Oy());

            Array3<float> zeta(oceanModelAdapter->Zeta().Nx(), oceanModelAdapter->Zeta().Ny(), oceanModelAdapter->Zeta().Nz(),
                               oceanModelAdapter->Zeta().Ox(), oceanModelAdapter->Zeta().Oy(), oceanModelAdapter->Zeta().Oz());
            Array4<float> u(
                    oceanModelAdapter->U().Nx(), oceanModelAdapter->U().Ny(), oceanModelAdapter->U().Nz(), oceanModelAdapter->U().N4(),
                    oceanModelAdapter->U().Ox(), oceanModelAdapter->U().Oy(), oceanModelAdapter->U().Oz(), oceanModelAdapter->U().O4());
            Array4<float> v(
                    oceanModelAdapter->V().Nx(), oceanModelAdapter->V().Ny(), oceanModelAdapter->V().Nz(), oceanModelAdapter->V().N4(),
                    oceanModelAdapter->V().Ox(), oceanModelAdapter->V().Oy(), oceanModelAdapter->V().Oz(), oceanModelAdapter->V().O4());
            Array4<float> w(
                    oceanModelAdapter->W().Nx(), oceanModelAdapter->W().Ny(), oceanModelAdapter->W().Nz(), oceanModelAdapter->W().N4(),
                    oceanModelAdapter->W().Ox(), oceanModelAdapter->W().Oy(), oceanModelAdapter->W().Oz(), oceanModelAdapter->W().O4());
            Array4<float> akt(
                    oceanModelAdapter->AKT().Nx(), oceanModelAdapter->AKT().Ny(), oceanModelAdapter->AKT().Nz(),
                    oceanModelAdapter->AKT().N4(),
                    oceanModelAdapter->AKT().Ox(), oceanModelAdapter->AKT().Oy(), oceanModelAdapter->AKT().Oz(),
                    oceanModelAdapter->AKT().O4());

            oceanTime.Load(oceanModelAdapter->OceanTime()());
            cout << "Wacomm::run oceanModelAdapter->OceanTime()("<< ocean_time_idx <<"): " << oceanModelAdapter->OceanTime()(ocean_time_idx) << endl;
            cout << "Wacomm::run oceanModelAdapter->Mask()(650,550):" << oceanModelAdapter->Mask()(650,550) << endl;
            cout << "Wacomm::run oceanModelAdapter->H()(650,550):" << oceanModelAdapter->H()(650,550) << endl;
            mask.Load(oceanModelAdapter->Mask()());
            cout << "Wacomm::run oceanModelAdapter->Mask()(650,550):" << oceanModelAdapter->Mask()(650,550) << endl;
            cout << "Wacomm::run oceanModelAdapter->H()(650,550):" << oceanModelAdapter->H()(650,550) << endl;
            exit(1);
            lonRad.Load(oceanModelAdapter->LonRad()());
            latRad.Load(oceanModelAdapter->LatRad()());
            depth.Load(oceanModelAdapter->Depth()());
            zeta.Load(oceanModelAdapter->Zeta()());
            u.Load(oceanModelAdapter->U()());
            v.Load(oceanModelAdapter->V()());
            w.Load(oceanModelAdapter->W()());
            akt.Load(oceanModelAdapter->AKT()());

            h.Load(oceanModelAdapter->H()());

            cout << "Wacomm::run oceanModelAdapter->H()(650,550):" << oceanModelAdapter->H()(650,550) << endl;
            cout << "Wacomm::run h(650,550):" << h(650,550) << endl;

            size_t nLocalParticles = pLocalParticles->size();

            size_t elementsPerThread = pLocalParticles->size() / ompMaxThreads;
            size_t sparePerThread = pLocalParticles->size() % ompMaxThreads;

            size_t thread_counts[ompMaxThreads];
            size_t thread_displs[ompMaxThreads];

            thread_counts[0] = (int)((elementsPerThread + sparePerThread));
            thread_displs[0] = 0;

            for (int tidx = 1; tidx < ompMaxThreads; tidx++) {
                thread_counts[tidx] = (int)(elementsPerThread);
                thread_displs[tidx]=thread_counts[0]+elementsPerThread*(tidx-1);
            }

            // Record start time
            auto startLocal = std::chrono::high_resolution_clock::now();

            #pragma omp parallel default(none) private(ompThreadNum) shared(thread_counts, thread_displs, pConfigData, pLocalParticles, ocean_time_idx, oceanTime, mask, lonRad,latRad, depth, h, zeta, u, v, w, akt)
            {

#ifdef USE_OMP
                ompThreadNum=omp_get_thread_num();
#endif

                struct config_data configData;
                memcpy(&configData,pConfigData,sizeof(configData));

                size_t first=thread_displs[ompThreadNum];
                size_t last=first+thread_counts[ompThreadNum];

                for (size_t idx = first; idx < last; idx++) {
                    pLocalParticles->at(idx).move(&configData, ocean_time_idx, oceanTime, mask, lonRad,
                                                  latRad, depth, h, zeta, u, v, w, akt);
                }

            }
            // Record end time
            auto finishLocal = std::chrono::high_resolution_clock::now();

#ifdef USE_MPI
            for (int idx=0;idx<localParticles.size();idx++) {
                recvbuf[idx]=localParticles.at(idx).data();
            }
            MPI_Gatherv(recvbuf.get(), send_counts.get()[world_rank],mpiParticleData,
                        sendbuf.get(),send_counts.get(),displs.get(),mpiParticleData,0,MPI_COMM_WORLD);
            MPI_Type_free(&mpiParticleData);
#endif

            if (world_rank==0) {
                std::chrono::duration<double> elapsedLocal = finishLocal - startLocal;
                double nParticlesPerSecondLocal = pLocalParticles->size() / elapsedLocal.count();
                LOG4CPLUS_INFO(logger, "Locally processed " << nParticles << " in " << elapsedLocal.count() << " seconds ("<< nParticlesPerSecondLocal <<" particles/second).");

#ifdef USE_MPI
                for (int idx=0;idx<particles->size();idx++) {
                    particles->at(idx).data(sendbuf[idx]);
                }
#endif
                // Remove dead particles
                particles->erase(std::remove_if(particles->begin(), particles->end(),
                                                [](const Particle &particle) { return !particle.isAlive(); }),
                                 particles->end());

                size_t nParticles = particles->size();
                LOG4CPLUS_INFO(logger, "Particles:" << nParticles);

                LOG4CPLUS_INFO(logger, "Evaluate concentration");

                #pragma omp parallel for default(none) shared(nParticles, ocean_time_idx, s_rho, eta_rho, xi_rho, conc)
                for (int idx=0; idx<nParticles;idx++) {
                    const Particle &particle = particles->at(idx);
                    if (particle.isAlive()) {
                        int k = (int) round(particle.K());
                        int j = (int) round(particle.J());
                        int i = (int) round(particle.I());
                        if (j>=0 && j<eta_rho && i>0 && i<xi_rho && k>=(-(int)s_rho+1) && k<=0) {
                            conc(ocean_time_idx, k, j, i) = conc(ocean_time_idx, k, j, i) + 1;
                        }
                    }
                }

                #pragma omp parallel for collapse(3) default(none) shared(ocean_time_idx, s_rho, eta_rho, xi_rho, conc)
                for (int k = -(int) s_rho + 1; k <= 0; k++) {
                    for (int j = 0; j < eta_rho; j++) {
                        for (int i = 0; i < xi_rho; i++) {
                            if (oceanModelAdapter->Mask()(j, i) != 1) {
                                conc(ocean_time_idx, k, j, i) = 1e37;
                            }
                        }
                    }
                }
            }
            auto finish = std::chrono::high_resolution_clock::now();

            if (world_rank==0) {
                std::chrono::duration<double> elapsed = finish - start;
                double nParticlesPerSecond = nParticles / elapsed.count();
                LOG4CPLUS_INFO(logger, "Processed " << nParticles << " in " << elapsed.count() << " seconds ("<< nParticlesPerSecond <<" particles/second).");
            }
        }



        if (world_rank==0) {
            Calendar calFinal;
            // Calculate the oceanTime at the end of calculations
            double finalOceanTime=oceanModelAdapter->OceanTime()[ocean_time-1]+config->Deltat();
            JulianDate::fromModJulian(finalOceanTime/86400, calFinal);

            if (!config->SaveHistory().empty()) {

                string historyFilename = config->HistoryRoot() + calFinal.asNCEPdate() ;
                LOG4CPLUS_INFO(logger, "Saving restart:" << historyFilename);
                if (config->SaveHistory() == "text") {
                    particles->saveAsTxt(historyFilename+ ".txt");
                } else if (config->SaveHistory() == "json") {
                    particles->saveAsJson(historyFilename+ ".json", finalOceanTime, oceanModelAdapter);
                } else if (config->SaveHistory() == "nc") {
                    particles->saveAsNetCDF(historyFilename + ".nc", finalOceanTime, oceanModelAdapter);
                }
            }

            string ncOutputFilename=config->NcOutputRoot()+cal.asNCEPdate()+".nc";
            LOG4CPLUS_INFO(logger, "Saving history:" << ncOutputFilename);
            if (config->EmbeddedHistory()) {
                particles->saveAsNetCDF(ncOutputFilename, finalOceanTime, oceanModelAdapter);
            }
            save(ncOutputFilename,conc, config->EmbeddedHistory());
        }
    }
};
Wacomm::~Wacomm() = default;

void Wacomm::save(string &fileName, Array4<float> &conc, bool add) {

    size_t ocean_time=oceanModelAdapter->OceanTime().Nx();
    size_t s_rho=oceanModelAdapter->SRho().Nx();
    size_t eta_rho=oceanModelAdapter->Mask().Nx();
    size_t xi_rho=oceanModelAdapter->Mask().Ny();

    NcFile::FileMode fileMode=NcFile::replace;
    if (add == true) {
        fileMode=NcFile::write;
    }



    LOG4CPLUS_INFO(logger,"Saving in: " << fileName);

    // Open the file for read access
    netCDF::NcFile dataFile(fileName, fileMode,NcFile::nc4);

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

