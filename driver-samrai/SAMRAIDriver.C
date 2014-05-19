//
// $Id: test_gradient.C 1832 2005-08-23 21:10:00Z bphilip $
// $Revision: 1832 $
// $Date: 2005-08-23 15:10:00 -0600 (Tue, 23 Aug 2005) $
//

#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <sys/stat.h>

// Application header.
#include "SAMRAIDriver.h"

// SAMRAI headers
#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/MemoryDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/appu/VisItDataWriter.h"
//#include "SAMRAI/solv/PETSc_SAMRAIVectorReal.h"

// SAMRUTILS headers
#include "testutils/SAMRBuilder.h"
#include "utilities/ProfilerApp.h"
#include "utilities/Utilities.h"

// SAMRSOLVERS headers
#include "time_integrators/TimeIntegratorParameters.h"
#include "time_integrators/TimeIntegratorFactory.h"

// Local headers
#include "pixie3dApplication.h"
#include "pixie3dApplicationParameters.h"

extern "C"{
#include "assert.h"
}

/*
 * Ghost cell width for variables that need them.
 */

#define GHOST_CELL_WIDTH (1)


int main( int argc, char *argv[] ) 
{
  SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);
  SAMRAI::tbox::SAMRAIManager::initialize();
  SAMRAI::tbox::SAMRAIManager::startup();
  const SAMRAI::tbox::SAMRAI_MPI& mpi(SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld());

  //PetscInitializeNoArguments();
  //PetscInitializeFortran();

  Utilities::setAllErrorHandlers( );

  // This extra code block is used to scope some temporaries that are
  // created, it forces the destruction before the manager is shutdown.
  {

    PROFILE_START("MAIN");
    std::string timer_results = "pixie3d.samrai";
    std::string input_file;
    std::string log_file;
   
    int rank = mpi.getRank();
    int rank2;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank2);
    TBOX_ASSERT(rank==rank2);

    // Process command line arguments and dump to log file.
    SAMRUtils::SAMRBuilder::processCommandLine(argc, argv, input_file, log_file);

    // Create the log file
    SAMRAI::tbox::PIO::logAllNodes(log_file);

    // Create input database and parse all data in input file.
    boost::shared_ptr<SAMRAI::tbox::MemoryDatabase> input_db(new SAMRAI::tbox::MemoryDatabase("input_db"));
    SAMRAI::tbox::InputManager::getManager()->parseInputFile(input_file, input_db);
    boost::shared_ptr<SAMRAI::tbox::Database> main_db = input_db->getDatabase("Main");
    boost::shared_ptr<SAMRAI::tbox::Database>  tag_db = input_db->getDatabase("StandardTagAndInitialize");

    // Get the output parameters
    bool use_visit = main_db->getBool("use_visit");
    double dt_save = 1e100;             // Default maximum time between saves
    int plot_interval = 0x7FFFFFFF;     // Default maximum number of iterators between saves
    std::vector<double> save_times;     // Default times to force a save
    std::string write_path = "output";  // Default path to save the data
    long int max_saves = 10000;         // Default maximum number of saves
    if ( use_visit ) {
        if ( main_db->keyExists("output_path") )
            write_path = main_db->getString("output_path");
        if ( main_db->keyExists("plot_interval") )
            plot_interval = main_db->getInteger("plot_interval");
        if ( main_db->keyExists("save_times") )
            save_times = main_db->getDoubleVector("save_times");
        if ( main_db->keyExists("dt_save") )
            dt_save = main_db->getDouble("dt_save");
        if ( main_db->keyExists("max_saves") )
            max_saves = (long int) main_db->getDouble("max_saves");
        timer_results = write_path + "/" + timer_results;
    }
    int save_debug = 0;
    if ( main_db->keyExists("save_debug") ) 
        save_debug = main_db->getInteger("save_debug");
    std::string debug_name = "debugFile";
    if ( main_db->keyExists("debug_name") )
        debug_name = main_db->getString("debug_name");

    // Options for regridding
    int regrid_interval = 0;
    if (main_db->keyExists("regrid_interval"))
        regrid_interval = main_db->getInteger("regrid_interval");

    // Create an empty pixie3dApplication (needed to create the StandardTagAndInitialize)
    const SAMRAI::tbox::Dimension dim(3);
    boost::shared_ptr<SAMRAI::pixie3dApplication> application( new SAMRAI::pixie3dApplication() );

    // Create the patch hierarchy
    boost::shared_ptr<SAMRAI::mesh::StandardTagAndInitStrategy> object = application;
    boost::shared_ptr<SAMRAI::mesh::GriddingAlgorithm> gridding_algorithm;
    boost::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy = SAMRUtils::SAMRBuilder::buildHierarchy(input_db,object,gridding_algorithm);
    TBOX_ASSERT(dim==hierarchy->getDim());

    // Check the initial hierarchy to see the patch distribution
    for ( int ln=0; ln<hierarchy->getNumberOfLevels(); ln++ ) {
        boost::shared_ptr<SAMRAI::hier::PatchLevel> level = hierarchy->getPatchLevel(ln);
        int N_local = level->getLocalNumberOfPatches();
        int N_global = level->getGlobalNumberOfPatches();
        SAMRAI::tbox::plog << ln << ": " << N_local << " of " << N_global << std::endl;
    }

    // Initialize the writer
    boost::shared_ptr<SAMRAI::appu::VisItDataWriter> visit_writer(
        new SAMRAI::appu::VisItDataWriter(dim,"pixie3d visualizer",write_path) );


    // Initialize the application
    boost::shared_ptr<SAMRAI::tbox::Database > pixie3d_db = input_db->getDatabase("pixie3d");
    boost::shared_ptr<SAMRAI::pixie3dApplicationParameters> application_parameters(
        new SAMRAI::pixie3dApplicationParameters(pixie3d_db) );
    application_parameters->d_hierarchy = hierarchy;
    application_parameters->d_VizWriter = visit_writer;
    application->initialize( application_parameters );

    // Initialize x0
    double t0 = 0.0;
    application->setInitialConditions(t0);

    // Perform a regrid to make sure the tagging is all set correctly
    std::vector<int> tag_buffer(20,2);
    if ( regrid_interval > 0 ) {
        boost::shared_ptr<mesh::StandardTagAndInitialize> error_detector = 
            boost::dynamic_pointer_cast<mesh::StandardTagAndInitialize>(gridding_algorithm->getTagAndInitializeStrategy());
        TBOX_ASSERT(error_detector!=NULL);
        error_detector->turnOnGradientDetector(t0);
        error_detector->turnOffRefineBoxes(t0);
        tag_buffer = std::vector<int>(20,2);    // GradientDetector or RichardsonExtrapolation is used, use a default tag buffer of 2
        for (int ln = 0; hierarchy->levelCanBeRefined(ln); ln++) {
            tbox::pout << "Regridding from level " << ln << std::endl;
            gridding_algorithm->regridAllFinerLevels(ln,tag_buffer,t0,0);
            application->setInitialConditions(t0);
        }
    }

    // initialize a time integrator parameters object
    boost::shared_ptr<tbox::Database> ti_db = input_db->getDatabase("TimeIntegrator");
    boost::shared_ptr<SAMRSolvers::TimeIntegratorParameters> timeIntegratorParameters(
       new SAMRSolvers::TimeIntegratorParameters(ti_db) );
    timeIntegratorParameters->d_operator = application;
    timeIntegratorParameters->d_ic_vector = application->get_ic();
    timeIntegratorParameters->d_vizWriter = visit_writer;
    double dt_min = ti_db->getDoubleWithDefault("min_dt",0.0);
    double dt_max = ti_db->getDoubleWithDefault("max_dt",1e100);
   
    // create a time integrator
    boost::shared_ptr<SAMRSolvers::TimeIntegratorFactory> tiFactory( new SAMRSolvers::TimeIntegratorFactory() );
    boost::shared_ptr<SAMRSolvers::TimeIntegrator> timeIntegrator = tiFactory->createTimeIntegrator(timeIntegratorParameters);

    // Create x(t)
    boost::shared_ptr< solv::SAMRAIVectorReal<double> > x_t;
    x_t = timeIntegrator->getCurrentSolution();

    // Write the data
    if ( use_visit )
        visit_writer->writePlotData(hierarchy,0,0);
    PROFILE_SAVE(timer_results);
    FILE *debug_file = NULL;
    if ( save_debug>0 ) {
        if ( rank==0 )
            debug_file = fopen( debug_name.c_str() , "wb" );
        application->writeDebugData(debug_file,0,0,save_debug);
    }

    // Loop through time
    bool first_step = true;
    double current_time = 0.0;
    double dt = 0.01;
    double final_time = timeIntegrator->getFinalTime();
    long int timestep = 0;
    double last_save_time = current_time;
    long int last_save_it = timestep;
    int it_save_time = 0;
    while ( current_time < final_time ) {
        current_time = timeIntegrator->getCurrentTime();
        
        // try and take a step
        PROFILE_START("advanceSolution");
        timeIntegrator->advanceSolution(dt, first_step);
        PROFILE_STOP("advanceSolution");
        // check if the computed approximation is acceptable for the timestep taken
        bool solnAcceptable = timeIntegrator->checkNewSolution();

        if(solnAcceptable) {
            first_step = false;
            timeIntegrator->updateSolution();
            current_time = timeIntegrator->getCurrentTime();
            tbox::pout << "Advanced solution to time : " << current_time << std::endl;
            timestep++;

            // Determine if we want to save the data
            bool save_now = false;
            if ( timestep-last_save_it >= plot_interval )
                save_now = true;
            if ( current_time-last_save_time >= dt_save )
                save_now = true;
            if ( (int) save_times.size() > it_save_time ) {
                if ( current_time >= save_times[it_save_time] ) {
                    save_now = true;
                    it_save_time++;
                }
            }
            if ( save_now ) {
                if ( use_visit )
                    visit_writer->writePlotData(hierarchy, timestep, current_time);
                if ( save_debug>0 )
                    application->writeDebugData(debug_file,timestep,current_time,save_debug);
                PROFILE_SAVE(timer_results);
                last_save_it = timestep;
                last_save_time = current_time;
            }

            // If desired, regrid patch hierarchy and reset vector weights.
            if ( regrid_interval>0 && ((timestep+1)%regrid_interval)==0 ) {
                tbox::pout << " Regridding ..." << std::endl;
                application->registerVector( x_t );
                gridding_algorithm->regridAllFinerLevels( 0, tag_buffer, current_time, 0 );
                first_step = false;
                tbox::pout << "************************* Finished regrid *********************************" << std::endl;
            }

            //dt = timeIntegrator->getNextDt(solnAcceptable);
            dt = application->getExpdT();
            dt = std::max(dt,dt_min);
            dt = std::min(dt,dt_max);
        
            tbox::pout << "Estimating next time step : " << dt << std::endl;

        } else {
            tbox::pout << "Failed to advance solution past time : " << timeIntegrator->getCurrentTime() << 
                ", current time step: " << timeIntegrator->getCurrentDt() << ", recomputing timestep ..." << std::endl;
			TBOX_ERROR("Error advancing solution");
        }

    }
    if ( debug_file != NULL )
        fclose(debug_file);

    // Barrier to make sure all processors have finished
    mpi.Barrier();
    // Delete the time integrator
    tiFactory.reset();
    timeIntegratorParameters.reset();
    // Delete the application
    application.reset();
    application_parameters.reset();
    PROFILE_STOP("MAIN");
    PROFILE_SAVE(timer_results);
    
  } // End code block
  tbox::pout << "Finializing PETSc, MPI and SAMRAI" << std::endl;

  // That's all, folks!
  mpi.Barrier();
  //PetscFinalize();
  SAMRAI::tbox::SAMRAIManager::shutdown();
  SAMRAI::tbox::SAMRAIManager::finalize();
  SAMRAI::tbox::SAMRAI_MPI::finalize(); 
  return(0);
}



