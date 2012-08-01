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
#include "SAMRAI/algs/ImplicitIntegrator.h"
#include "SAMRAI/solv/SNES_SAMRAIContext.h"

// Local headers
#include "ImplicitPixie3dApplication.h"
#include "ImplicitPixie3dApplicationParameters.h"

extern "C"{
#include "assert.h"
}

#include "petsc.h"

/*
 * Ghost cell width for variables that need them.
 */

#define GHOST_CELL_WIDTH (1)


int main( int argc, char *argv[] ) 
{
  SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);
  SAMRAI::tbox::SAMRAIManager::initialize();
  SAMRAI::tbox::SAMRAIManager::startup();
  const int maxPatchDataEntries = 2000;
  SAMRAI::tbox::SAMRAIManager::setMaxNumberPatchDataEntries(maxPatchDataEntries);
  const SAMRAI::tbox::SAMRAI_MPI& mpi(SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld());
  Utilities::setMPIErrorHandler( mpi );

  int ierr = PetscInitialize(&argc, &argv, PETSC_NULL,PETSC_NULL);
  //PetscInitializeNoArguments();
  //PetscInitializeFortran();

  // This extra code block is used to scope some temporaries that are
  // created, it forces the destruction before the manager is shutdown.
  {
    PROFILE_START("MAIN");
    std::string timer_results = "pixie3d.samrai";
    std::string input_file;
    std::string log_file;
   
    int rank = mpi.getRank();

    // Process command line arguments and dump to log file.
    SAMRUtils::SAMRBuilder::processCommandLine(argc, argv, input_file, log_file);

    // Create the log file
    SAMRAI::tbox::PIO::logAllNodes(log_file);

    // Create input database and parse all data in input file.
    SAMRAI::tbox::Pointer<SAMRAI::tbox::MemoryDatabase> input_db(new SAMRAI::tbox::MemoryDatabase("input_db"));
    SAMRAI::tbox::InputManager::getManager()->parseInputFile(input_file, input_db);
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> main_db = input_db->getDatabase("Main");
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>  tag_db = input_db->getDatabase("StandardTagAndInitialize");

    // Get the output parameters
    bool use_visit = main_db->getBool("use_visit");
    double dt_save = 1e100;             // Default maximum time between saves
    int plot_interval = 0x7FFFFFFF;     // Default maximum number of iterators between saves
    tbox::Array<double> save_times;     // Default times to force a save
    std::string write_path = "output";  // Default path to save the data
    long int max_saves = 10000;         // Default maximum number of saves
    if ( use_visit ) {
        if ( main_db->keyExists("output_path") )
            write_path = main_db->getString("output_path");
        if ( main_db->keyExists("plot_interval") )
            plot_interval = main_db->getInteger("plot_interval");
        if ( main_db->keyExists("save_times") )
            save_times = main_db->getDoubleArray("save_times");
        if ( main_db->keyExists("dt_save") )
            dt_save = main_db->getDouble("dt_save");
        if ( main_db->keyExists("max_saves") )
            max_saves = (long int) main_db->getDouble("max_saves");
    }
    int save_debug = 0;
    if ( main_db->keyExists("save_debug") ) 
        save_debug = main_db->getInteger("save_debug");
    std::string debug_name = "debugFile";
    if ( main_db->keyExists("debug_name") )
        debug_name = main_db->getString("debug_name");
    timer_results = write_path + "/" + timer_results;
	
    // Options for printing residuals
    bool print_nonlinear_residuals = false;
    if (main_db->keyExists("print_nonlinear_residuals"))
        print_nonlinear_residuals = main_db->getBool("print_nonlinear_residuals");
    bool print_linear_residuals = false;
    if (main_db->keyExists("print_linear_residuals"))
        print_linear_residuals = main_db->getBool("print_linear_residuals");
    
    // Options for regridding
    int regrid_interval = 0;
    if (main_db->keyExists("regrid_interval"))
        regrid_interval = main_db->getInteger("regrid_interval");

    // Create an empty pixie3dApplication (needed to create the StandardTagAndInitialize)
    const SAMRAI::tbox::Dimension dim(3);
    tbox::Pointer<SAMRAI::Pixie3d::ImplicitPixie3dApplication> application( new SAMRAI::Pixie3d::ImplicitPixie3dApplication() );

    // Create the patch hierarchy
    SAMRAI::tbox::Pointer<SAMRAI::mesh::StandardTagAndInitStrategy> object = application;
    SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm> gridding_algorithm;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy> hierarchy = SAMRUtils::SAMRBuilder::buildHierarchy(input_db,object,gridding_algorithm);
    TBOX_ASSERT(dim==hierarchy->getDim());

    // Check the initial hierarchy to see the patch distribution
    for ( int ln=0; ln<hierarchy->getNumberOfLevels(); ln++ ) {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel> level = hierarchy->getPatchLevel(ln);
        int N_local = level->getLocalNumberOfPatches();
        int N_global = level->getGlobalNumberOfPatches();
        SAMRAI::tbox::plog << ln << ": " << N_local << " of " << N_global << std::endl;
    }

    // Initialize the writer
    appu::VisItDataWriter* visit_writer;
    visit_writer = new appu::VisItDataWriter( dim, "pixie3d visualizer", write_path );


    // Initialize the application
    tbox::Pointer< tbox::Database > pixie3d_db = input_db->getDatabase("pixie3d");
    SAMRAI::Pixie3d::ImplicitPixie3dApplicationParameters* application_parameters = new SAMRAI::Pixie3d::ImplicitPixie3dApplicationParameters(pixie3d_db);
    application_parameters->d_hierarchy = hierarchy;
    application_parameters->d_VizWriter = visit_writer;
    application->initialize( application_parameters );

    // Initialize x0
    double t0 = 0.0;
    application->setInitialConditions(t0);

    // Perform a regrid to make sure the tagging is all set correctly
    tbox::Array<int> tag_buffer(20,2);
    if ( regrid_interval > 0 ) {
        tbox::Pointer<mesh::StandardTagAndInitialize> error_detector = gridding_algorithm->getTagAndInitializeStrategy();
        TBOX_ASSERT(!error_detector.isNull());
        error_detector->turnOnGradientDetector();
        error_detector->turnOffRefineBoxes();
        tag_buffer = tbox::Array<int>(20,2);    // GradientDetector or RichardsonExtrapolation is used, use a default tag buffer of 2
        for (int ln = 0; hierarchy->levelCanBeRefined(ln); ln++) {
            tbox::pout << "Regridding from level " << ln << std::endl;
            gridding_algorithm->regridAllFinerLevels(ln,t0,tag_buffer);
            application->setInitialConditions(t0);
        }
    }

    //application->createPreconditioner();
    
    solv::SNES_SAMRAIContext* snes_solver = NULL;
    
    snes_solver = new solv::SNES_SAMRAIContext("SNESSolver",
                           input_db->getDatabase("SNESSolver"),
                           application);

    tbox::pout << "Created nonlinear solver " << std::endl;
    
    input_db->getDatabase("SNESSolver")->printClassData(tbox::plog);

    algs::ImplicitIntegrator* timeIntegrator = 
        new algs::ImplicitIntegrator("ImplicitIntegrator",
                   input_db->getDatabase("ImplicitIntegrator"),
                   application,
                   snes_solver,
                   hierarchy);

    tbox::pout << "Created implicit time integrator " << std::endl;
    input_db->getDatabase("ImplicitIntegrator")->printClassData(tbox::plog);

    /*
     * Initialize implicit integrator.  Then, loop over a sequence of time 
     * steps until the final simulation time is reached, or the maximum 
     * number of steps is exceeded.
     */
    timeIntegrator->initialize();
    
    tbox::pout << "Initialized implicit time integrator " << std::endl;

    std::string gmresOrthogonalizationMethod="modifiedgramschmidt";
    
    snes_solver->setGMRESOrthogonalizationMethod(gmresOrthogonalizationMethod);
    
#if PETSC_VERSION_(3,0,0)   
    ierr = KSPSetPreconditionerSide(snes_solver->getKrylovSolver(),
                    PC_RIGHT);
    
    PETSC_SAMRAI_ERROR(ierr);
#else
    ierr = KSPSetPCSide(snes_solver->getKrylovSolver(),
            PC_RIGHT);
    
    PETSC_SAMRAI_ERROR(ierr);
#endif
    
   
    if ( print_nonlinear_residuals ) {
        ierr = SNESMonitorSet(snes_solver->getSNESSolver(),
                  SNESMonitorDefault,
                  PETSC_NULL,
                  PETSC_NULL);  CHKERRQ(ierr);
    }
    
    if ( print_linear_residuals ) {
        ierr = KSPMonitorSet(snes_solver->getKrylovSolver(),
                 KSPMonitorDefault, 
                 PETSC_NULL, 
                 PETSC_NULL);  CHKERRQ(ierr);
    }

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
    double dt = timeIntegrator->getCurrentDt();
    dt = timeIntegrator->getNextDt(true,0);
    double final_time = timeIntegrator->getFinalTime();
    long int timestep = 0;
    double last_save_time = current_time;
    long int last_save_it = timestep;
    int it_save_time = 0;
    // the next two variables are to keep statistics on the
    // total number of linear and nonlinear solves
    int total_nonlinear_itns = 0;
    int total_linear_itns = 0;
    while ( current_time < final_time ) {
        current_time = timeIntegrator->getCurrentTime();
        tbox::pout << "Timestep " << timestep  << ", current time: " << current_time << std::endl;
        
        // try and take a step
        PROFILE_START("advanceSolution");
        int solver_retcode = timeIntegrator->advanceSolution(dt, first_step);
        PROFILE_STOP("advanceSolution");
        int nonlinear_itns, linear_itns;
        nonlinear_itns = snes_solver->getNumberOfNonlinearIterations();
        linear_itns = snes_solver->getTotalNumberOfLinearIterations();
        tbox::pout << "Completed nonlinear solve" << std::endl;
        total_nonlinear_itns += nonlinear_itns;
        total_linear_itns += linear_itns;

        // check if the computed approximation is acceptable for the timestep taken
        bool solnAcceptable = timeIntegrator->checkNewSolution(solver_retcode);

        if(solnAcceptable) {
            first_step = false;

            // If desired, regrid patch hierarchy and reset vector weights.
            if ( regrid_interval>0 && ((timestep+1)%regrid_interval)==0 ) {

                tbox::pout << " Regridding ..." << std::endl;
                gridding_algorithm->regridAllFinerLevels( 0, current_time, tag_buffer );
                tbox::pout << "************************* Finished regrid *********************************" << std::endl;
   
                snes_solver->resetSolver(0, hierarchy->getFinestLevelNumber());
                
                snes_solver->setGMRESOrthogonalizationMethod(gmresOrthogonalizationMethod);
        
                #if PETSC_VERSION_(3,0,0)   
                    ierr = KSPSetPreconditionerSide( snes_solver->getKrylovSolver(), PC_RIGHT );
                    PETSC_SAMRAI_ERROR(ierr);
                #else
                    ierr = KSPSetPCSide( snes_solver->getKrylovSolver(), PC_RIGHT );
                    PETSC_SAMRAI_ERROR(ierr);
                #endif
        
                if ( print_nonlinear_residuals ) {
                    ierr = SNESMonitorSet(snes_solver->getSNESSolver(),
                        SNESMonitorDefault,
                        PETSC_NULL,
                        PETSC_NULL);  CHKERRQ(ierr);
                }
        
                if ( print_linear_residuals ) {
                    ierr = KSPMonitorSet(snes_solver->getKrylovSolver(),
                        KSPMonitorDefault, 
                        PETSC_NULL, 
                        PETSC_NULL);  CHKERRQ(ierr);
                }
                solver_retcode = timeIntegrator->advanceSolution(dt, first_step);

            }
			
            current_time = timeIntegrator->updateSolution();
            tbox::pout << "Advanced solution to time : " << current_time << std::endl;
            timestep++;

            // Determine if we want to save the data
            bool save_now = false;
            if ( timestep-last_save_it >= plot_interval )
                save_now = true;
            if ( current_time-last_save_time >= dt_save )
                save_now = true;
            if ( save_times.size() > it_save_time ) {
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

            dt = timeIntegrator->getNextDt(solnAcceptable, solver_retcode);
            tbox::pout << "\n\n++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
            tbox::pout << "At end of timestep # " << timestep - 1 << std::endl;
            tbox::pout << " Iterations:  nonlinear " << nonlinear_itns << std::endl;
            tbox::pout << "              linear:   " << linear_itns << std::endl;
            tbox::pout << "Simulation time is " << current_time << std::endl;
            tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++\n\n" << std::endl;
            tbox::pout << "Estimating next time step : " << dt << std::endl;
        } else {
            tbox::pout << "Failed to advance solution past time : " << timeIntegrator->getCurrentTime() << ", current time step: " << timeIntegrator->getCurrentDt() << ", recomputing timestep ..." << std::endl;
			TBOX_ERROR("Error advancing solution");
        }

    }
    if ( debug_file != NULL )
        fclose(debug_file);

    // Barrier to make sure all processors have finished
    mpi.Barrier();
    // Delete the application
    application.setNull();
    delete application_parameters;
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



