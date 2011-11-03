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

  // This extra code block is used to scope some temporaries that are
  // created, it forces the destruction before the manager is shutdown.
  {

    std::string input_file;
    std::string log_file;
    int plot_interval=0;
    std::string write_path;
   
    int rank = mpi.getRank();

    // Process command line arguments and dump to log file.
    SAMRAI::SAMRBuilder::processCommandLine(argc, argv, input_file, log_file);

    // Create the log file
    SAMRAI::tbox::PIO::logOnlyNodeZero(log_file);

    // Create input database and parse all data in input file.
    SAMRAI::tbox::Pointer<SAMRAI::tbox::MemoryDatabase> input_db(new SAMRAI::tbox::MemoryDatabase("input_db"));
    SAMRAI::tbox::InputManager::getManager()->parseInputFile(input_file, input_db);
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> main_db = input_db->getDatabase("Main");
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>  tag_db = input_db->getDatabase("StandardTagAndInitialize");
    double dt_save;
    if ( main_db->keyExists("output_path") ) {
        write_path = main_db->getString("output_path");
    } else {
        write_path = "output";
    }
    if ( main_db->keyExists("dt_save") )
        dt_save = main_db->getDouble("dt_save");
    else
        dt_save = 1.0;
    int save_debug = 0;
    if ( main_db->keyExists("save_debug") ) 
        save_debug = main_db->getInteger("save_debug");
    std::string debug_name;
    if ( main_db->keyExists("debug_name") )
        debug_name = main_db->getString("debug_name");
    else
        debug_name = "debugFile";

    // Create an empty pixie3dApplication (needed to create the StandardTagAndInitialize)
    const SAMRAI::tbox::Dimension dim(3);
    tbox::Pointer<SAMRAI::pixie3dApplication> application( new SAMRAI::pixie3dApplication() );

    // Create the patch hierarchy
    SAMRAI::tbox::Pointer<SAMRAI::mesh::StandardTagAndInitStrategy> object = application;
    SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm> gridding_algorithm;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy> hierarchy = SAMRAI::SAMRBuilder::buildHierarchy(input_db,object,gridding_algorithm);
    TBOX_ASSERT(dim==hierarchy->getDim());

    // Check that the initial hierarchy has at least one patch per processor
    int N_patches = 0;
    for ( int ln=0; ln<hierarchy->getNumberOfLevels(); ln++ ) {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel> level = hierarchy->getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel::Iterator p(level); p; p++)
            N_patches++;
    }
    if ( N_patches==0 )
        TBOX_ERROR("One or more processors does not have any patches");

    // Initialize the writer
    appu::VisItDataWriter* visit_writer;
    visit_writer = new appu::VisItDataWriter( dim, "pixie3d visualizer", write_path );


    // Initialize the application
    tbox::Pointer< tbox::Database > pixie3d_db = input_db->getDatabase("pixie3d");
    SAMRAI::pixie3dApplicationParameters* application_parameters = new SAMRAI::pixie3dApplicationParameters(pixie3d_db);
    application_parameters->d_hierarchy = hierarchy;
    application_parameters->d_VizWriter = visit_writer;
    application->initialize( application_parameters );

    // Initialize x0
    double t0 = 0.0;
    application->setInitialConditions(t0);

    // Perform a regrid to make sure the tagging is all set correctly
    tbox::Pointer<mesh::StandardTagAndInitialize> error_detector = gridding_algorithm->getTagAndInitializeStrategy();
    TBOX_ASSERT(!error_detector.isNull());
    error_detector->turnOnGradientDetector();
    error_detector->turnOffRefineBoxes();
    tbox::Array<int> tag_buffer;
    if ( error_detector->refineUserBoxInputOnly() )
        tag_buffer = tbox::Array<int>(20,0);    // Only user-defined refinement boxes are used, no tag buffer is necessary
    else
        tag_buffer = tbox::Array<int>(20,2);    // GradientDetector or RichardsonExtrapolation is used, use a default tag buffer of 2
    gridding_algorithm->regridAllFinerLevels(0,t0,tag_buffer);
    application->setInitialConditions(t0);

    // initialize a time integrator parameters object
    tbox::Pointer<tbox::Database> ti_db = input_db->getDatabase("TimeIntegrator");
    SAMRSolvers::TimeIntegratorParameters *timeIntegratorParameters = new SAMRSolvers::TimeIntegratorParameters(ti_db);
    timeIntegratorParameters->d_operator = application;
    timeIntegratorParameters->d_ic_vector = application->get_ic();
    timeIntegratorParameters->d_vizWriter = visit_writer;
   
    // create a time integrator
    SAMRSolvers::TimeIntegratorFactory *tiFactory = new SAMRSolvers::TimeIntegratorFactory();
    std::auto_ptr<SAMRSolvers::TimeIntegrator> timeIntegrator = tiFactory->createTimeIntegrator(timeIntegratorParameters);

    // Create x(t)
    tbox::Pointer< solv::SAMRAIVectorReal<double> > x_t;
    x_t = timeIntegrator->getCurrentSolution();

    tbox::Pointer<tbox::Database> plot_db = input_db->getDatabase("Plotting");

    if (plot_db->keyExists("plot_interval")) {
        plot_interval = plot_db->getInteger("plot_interval");
    }

    // Write the data
    visit_writer->writePlotData(hierarchy,0,0);
    FILE *debug_file = NULL;
    if ( save_debug>0 ) {
        if ( rank==0 )
            debug_file = fopen( debug_name.c_str() , "wb" );
        application->writeDebugData(debug_file,0,0,save_debug);
    }

    // Loop through time
    bool first_step = true;
    double dt = 1.0;

    double current_time = 0.0;
    double final_time = timeIntegrator->getFinalTime();
   
    int iteration_num = 0;
    while (current_time < final_time) {
        iteration_num++;
        current_time = timeIntegrator->getCurrentTime();
        timeIntegrator->advanceSolution(dt, first_step);
        bool solnAcceptable = timeIntegrator->checkNewSolution();

        if(solnAcceptable) {
            first_step = false;
            timeIntegrator->updateSolution();
            current_time = timeIntegrator->getCurrentTime();
            tbox::pout << "Advanced solution to time : " << current_time << std::endl;

            if ( (plot_interval > 0) && ((iteration_num % plot_interval) == 0) )  {
                if ( save_debug>0 )
                    application->writeDebugData(debug_file,iteration_num,current_time,save_debug);
                visit_writer->writePlotData(hierarchy, iteration_num, current_time);
            }
        } else {
            tbox::pout << "Failed to advance solution past time : " << timeIntegrator->getCurrentTime() << ", current time step: " << timeIntegrator->getCurrentDt() << ", recomputing timestep ..." << std::endl;
        }

        //dt = timeIntegrator->getNextDt(solnAcceptable);
        dt = application->getExpdT();
        
        tbox::pout << "Estimating next time step : " << dt << std::endl;
    }

    if ( debug_file != NULL )
        fclose(debug_file);

    // Barrier to make sure all processors have finished
    mpi.Barrier();
    // Delete the time integrator
    delete tiFactory;
    delete timeIntegratorParameters;
    // Delete the application
    application.setNull();
    delete application_parameters;

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



