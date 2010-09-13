//
// $Id: test_gradient.C 1832 2005-08-23 21:10:00Z bphilip $
// $Revision: 1832 $
// $Date: 2005-08-23 15:10:00 -0600 (Tue, 23 Aug 2005) $
//

#include "SAMRAI_config.h"

#include <iostream>

#include <fstream>
using namespace std;

#include <sys/stat.h>

/*
 * SAMRAI headers.
 */
#include "BoundaryBox.h"
#include "BoxArray.h"
#include "BoxList.h"
#include "BergerRigoutsos.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellVariable.h"
#include "FaceData.h"
#include "FaceVariable.h"
#include "GriddingAlgorithm.h"
#include "Index.h"
#include "tbox/InputDatabase.h"
#include "tbox/InputManager.h"
#include "IntVector.h"
#include "LoadBalancer.h"
#include "tbox/SAMRAI_MPI.h"
#include "Patch.h"
#include "PatchGeometry.h"
#include "tbox/PIO.h"
#include "RefineAlgorithm.h"
#include "RefineOperator.h"
#include "tbox/SAMRAIManager.h"
#include "StandardTagAndInitialize.h"
#include "tbox/Utilities.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "PETSc_SAMRAIVectorReal.h"

#include "pixie3dApplication.h"
#include "pixie3dApplicationParameters.h"
#include "TimeIntegratorParameters.h"
#include "TimeIntegratorFactory.h"
#include "VisItDataWriter.h"
#include <sstream>
#include <string>

extern "C"{
#include "assert.h"
}
/*
 * Application header.
 */
#include "BogusTagAndInitStrategy.h"
#include "SAMRAIDriver.h"

/*
 * Ghost cell width for variables that need them.
 */

#define GHOST_CELL_WIDTH (1)

#define USE_MARK

int main( int argc, char *argv[] ) 
{
   string input_file;
   string log_file;
   int plot_interval=0;
   string write_path;
   
   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::startup();
   int ierr = PetscInitializeNoArguments();
   PetscInitializeFortran();

   /*
    * Process command line arguments and dump to log file.
    */
   processCommandLine(argc, argv, input_file, log_file);

   tbox::PIO::logOnlyNodeZero(log_file);

   /*
    * Create input database and parse all data in input file.  This
    * parsing allows us to subsequently extract individual sections.
    */
   tbox::Pointer<tbox::Database> input_db = new tbox::InputDatabase("input_db");
   tbox::InputManager::getManager()->parseInputFile(input_file, input_db);
   tbox::Pointer< tbox::Database > main_db = input_db->getDatabase("Main");
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

   /*
    * Create an application object.  Some TagAndInitStrategy must be
    * provided in order to build an object that specifies cells that
    * need refinement.  Here an empty object is provided, since a
    * prescribed set of refinement regions are read in from the input
    * file; it would be a useful exercise to fill in the
    * applyGradientDetector method and to generate custom refinement
    * regions.
    */
   BogusTagAndInitStrategy* test_object = new BogusTagAndInitStrategy();

   tbox::Pointer<tbox::Database> griddingDb = input_db->getDatabase("GriddingAlgorithm");
   tbox::Pointer<tbox::Database> ratioDb = griddingDb->getDatabase("ratio_to_coarser");
   // Create dummy variable to allow larger ghost cell widths
   hier::IntVector<NDIM> ghost = hier::IntVector<NDIM>::IntVector(ratioDb->getIntegerArray("level_1"));
   //   ghost(2) = 1;
   tbox::Pointer< pdat::CellVariable<NDIM,double> > var = new pdat::CellVariable<NDIM,double>( "tmp", 1 );
   hier::VariableDatabase<NDIM>* var_db = hier::VariableDatabase<NDIM>::getDatabase();
   tbox::Pointer<hier::VariableContext> d_application_ctx = var_db->getContext("APPLICATION_SCRATCH");
   int var_id = var_db->registerVariableAndContext(var, d_application_ctx, ghost );

   // Create the AMR hierarchy and initialize it
   tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy;
   #ifdef USE_MARK
      // Get some of the databases
      tbox::Pointer< tbox::Database > grid_db = input_db->getDatabase("CartesianGeometry");
      tbox::Pointer< tbox::Database >  tag_db = input_db->getDatabase("StandardTagAndInitialize");
      tbox::Pointer< tbox::Database > load_db = input_db->getDatabase("LoadBalancer");
      tbox::Pointer< tbox::Database > grid_alg_db = input_db->getDatabase("GriddingAlgorithm");
      // Create the hierarchy
      tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry =
         new geom::CartesianGridGeometry<NDIM>("CartesianGeometry",grid_db);
      hierarchy = new hier::PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
      // Create the gridding algorithum
      tbox::Pointer<mesh::StandardTagAndInitialize<NDIM> > error_detector =
         new mesh::StandardTagAndInitialize<NDIM>( "CellTaggingMethod", test_object, tag_db );
      tbox::Pointer<mesh::BergerRigoutsos<NDIM> > box_generator = new mesh::BergerRigoutsos<NDIM>();
      tbox::Pointer<mesh::LoadBalancer<NDIM> > load_balancer =
         new mesh::LoadBalancer<NDIM>("UniformLoadBalance",load_db);
      tbox::Pointer<mesh::GriddingAlgorithm<NDIM> > gridding_algorithm =
         new mesh::GriddingAlgorithm<NDIM>("GriddingAlgorithm", grid_alg_db, error_detector, box_generator, load_balancer);
      // Create the coarsest level
      gridding_algorithm->makeCoarsestLevel(hierarchy, -1);
  #else
      initializeAMRHierarchy(input_db,test_object,hierarchy);
  #endif
 
   // Initialize the writer
   appu::VisItDataWriter<NDIM>* visit_writer;
   visit_writer = new appu::VisItDataWriter<NDIM>("pixie3d visualizer", write_path);


   // Create the application
   SAMRAI::pixie3dApplicationParameters* application_parameters = new SAMRAI::pixie3dApplicationParameters( input_db->getDatabase("pixie3d"));
   application_parameters->d_hierarchy = hierarchy;
   application_parameters->d_VizWriter = visit_writer;
   
   SAMRAI::pixie3dApplication* application  = new SAMRAI::pixie3dApplication( application_parameters );
   

   // Initialize x0
   application->setInitialConditions(0.0);

   // initialize a time integrator parameters object
   tbox::Pointer<tbox::Database> ti_db = input_db->getDatabase("TimeIntegrator");
   SAMRSolvers::TimeIntegratorParameters *timeIntegratorParameters = new SAMRSolvers::TimeIntegratorParameters(ti_db);
   timeIntegratorParameters->d_operator = application;
   timeIntegratorParameters->d_ic_vector = application->get_x();
   timeIntegratorParameters->d_vizWriter = visit_writer;
   
   // create a time integrator
   SAMRSolvers::TimeIntegratorFactory *tiFactory = new SAMRSolvers::TimeIntegratorFactory();
   std::auto_ptr<SAMRSolvers::TimeIntegrator> timeIntegrator = tiFactory->createTimeIntegrator(timeIntegratorParameters);

   // Create x(t)
   tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > x_t;
   x_t = timeIntegrator->getCurrentSolution();

   tbox::Pointer<tbox::Database> plot_db = input_db->getDatabase("Plotting");

   if (plot_db->keyExists("plot_interval")) 
     {
       plot_interval = plot_db->getInteger("plot_interval");

     }

   // Write the data
   visit_writer->writePlotData(hierarchy,0,0);

   // Loop through time
   bool first_step = true;
   double dt = 1.0;

   double current_time = 0.0;
   double final_time = timeIntegrator->getFinalTime();
   
   int iteration_num = 0;

   // this changes the set of applied boundary conditions
   // check with Luis if this is in the right place
   application->setBoundarySchedules(false);
   
   while (current_time < final_time)
     {       
       iteration_num++;
       current_time = timeIntegrator->getCurrentTime();
       timeIntegrator->advanceSolution(dt, first_step);
       bool solnAcceptable = timeIntegrator->checkNewSolution();

       if(solnAcceptable)
	 {
	   first_step = false;
	   timeIntegrator->updateSolution();
	   current_time = timeIntegrator->getCurrentTime();
	   tbox::pout << "Advanced solution to time : " << current_time << std::endl;

	   if ( (plot_interval > 0) && ((iteration_num % plot_interval) == 0) ) 
	     {
	       visit_writer->writePlotData(hierarchy, iteration_num, current_time);
	     }

	 }
       else
	 {
	   tbox::pout << "Failed to advance solution past time : " << timeIntegrator->getCurrentTime() << ", current time step: " << timeIntegrator->getCurrentDt() << ", recomputing timestep ..." << std::endl;
	 }

       dt = timeIntegrator->getNextDt(solnAcceptable);
       tbox::pout << "Estimating next time step : " << dt << std::endl;
     }



   // Delete the time integrator
   delete tiFactory;
   delete timeIntegratorParameters;
   // Delete the application
   delete application;
   delete application_parameters;
   // Delete the writer
   delete visit_writer;
   // Delete the hierarchy
   #ifdef USE_MARK
      //delete box_generator;
      //delete gridding_algorithm;
   #endif
   hierarchy.setNull();
   delete grid_geometry;
   // Delete the application object
   delete test_object;
   // Delete the input database
   main_db.setNull();
   input_db.setNull();
   #ifdef USE_MARK
      grid_db.setNull();
      tag_db.setNull();
      load_db.setNull();
      grid_alg_db.setNull();
   #endif
   tbox::InputManager::freeManager();
   delete input_db;
   tbox::pout << "Input database deleted\n";
   // Finalize PETsc, MPI, and SAMRAI
   tbox::SAMRAI_MPI::barrier();
   PetscFinalize();
   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAI_MPI::finalize();
   return(0);
}

/*
************************************************************************
*                                                                      *
*  Parse command line arguments, returning name of input file and log  *
*  file.                                                               *
*                                                                      *
************************************************************************
*/
void processCommandLine(int argc, 
                        char *argv[], 
                        string& input_file, 
                        string& log_file)
{
  if ( (argc != 3) ) {
    tbox::pout << "USAGE:  " << argv[0] << " <input file> <log file> " << endl;
    exit(-1);
  } else {
    input_file = argv[1];
    log_file = argv[2];
  }

  return;
}

/*
************************************************************************
*                                                                      *
* Generate a patch hierarchy from information specified in input       *
* database.                                                            *
*                                                                      *
************************************************************************
*/
void initializeAMRHierarchy(tbox::Pointer<tbox::Database> &input_db,
			    mesh::StandardTagAndInitStrategy<NDIM>* user_tagging_strategy,
			    tbox::Pointer<hier::PatchHierarchy<NDIM> > &hierarchy)
{
   /*
    * Create geometry object.  This specifies the index space of the
    * coarsest level, as well as its physical (Cartesian) coordinates.
    */
   tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry = 
      new geom::CartesianGridGeometry<NDIM>("CartesianGeometry",
                                input_db->getDatabase("CartesianGeometry"));

   /*
    * Create patch hierarchy.
    */
   hierarchy = new hier::PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);

   /* 
    * A mesh::GriddingAlgorithm<NDIM> is used to build the initial grid hierarchy.
    * Classes for tagging cells that need refinement, generation of
    * boxes from these tagged cells, and load balancing the grid
    * hierarchy are needed to build the mesh::GriddingAlgorithm<NDIM>.
    *
    * First build the object used to tag cells that need refinement.
    */
    tbox::Pointer<mesh::StandardTagAndInitialize<NDIM> > error_detector = 
	new mesh::StandardTagAndInitialize<NDIM>( 
	    "CellTaggingMethod", 
	    user_tagging_strategy, 
	    input_db->getDatabase("StandardTagAndInitialize"));
    
   /*
    * Next, specify the built-in Berger-Rigoutsos method for
    * generating boxes from the tagged cells.
    */
   tbox::Pointer<mesh::BergerRigoutsos<NDIM> > box_generator = new mesh::BergerRigoutsos<NDIM>(); 

   /*
    * Next, specify the built-in uniform load balancer to distribute
    * patches across processors.
    */
   tbox::Pointer<mesh::LoadBalancer<NDIM> > load_balancer =
      new mesh::LoadBalancer<NDIM>(input_db->getDatabase("LoadBalancer"));

   /*
    * Finally, build the grid generator, registering the above
    * strategies for tagging cells, generating boxes, and load
    * balancing the calculation.
    */
   tbox::Pointer<mesh::GriddingAlgorithm<NDIM> > gridding_algorithm =
      new mesh::GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                            input_db->getDatabase("GriddingAlgorithm"),
                            error_detector,
                            box_generator,
                            load_balancer);

   /*
    * Build an initial grid hierarchy.  Note that in this simple
    * example we do not buffer the refinement regions.
    */
   gridding_algorithm->makeCoarsestLevel(hierarchy, 0.0);
   
   bool done = false;
   bool initial_time = true;
   for (int ln = 0;
        gridding_algorithm->levelCanBeRefined(ln) && !done; 
        ln++) {
       gridding_algorithm->makeFinerLevel(hierarchy,
                                          0.0,
                                          initial_time,
                                          0);
       done = !(hierarchy->finerLevelExists(ln));
   }
}

