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
#include "tbox/MPI.h"
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
#include "ForwardEulerTimeIntegrator.h"
//#include "ExplicitPCTimeIntegrator.h"
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
/*#include "pims_local_struct.h"*/
/*
 * Ghost cell width for variables that need them.
 */

#define GHOST_CELL_WIDTH (1)

int main( int argc, char *argv[] ) 
{
   double max_error = 0.0;
   double l2_error = 0.0;
   string input_file;
   string log_file;
   
   tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy;

   tbox::MPI::init(&argc, &argv);
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
   string write_path;
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

   // Create dummy variable to allow larger goast cell widths
   hier::IntVector<NDIM> ghost = hier::IntVector<NDIM>::IntVector(2);
   ghost(2) = 1;
   tbox::Pointer< pdat::CellVariable<NDIM,double> > var = new pdat::CellVariable<NDIM,double>( "tmp", 1 );
   hier::VariableDatabase<NDIM>* var_db = hier::VariableDatabase<NDIM>::getDatabase();
   tbox::Pointer<hier::VariableContext> d_application_ctx = var_db->getContext("APPLICATION_SCRATCH");
   int var_id = var_db->registerVariableAndContext(var, d_application_ctx, ghost );

   // Create the AMR hierarchy and initialize it
   initializeAMRHierarchy(input_db,test_object,hierarchy);
 
   // Create the application
   pixie3dApplicationParameters* application_parameters = new pixie3dApplicationParameters();
   application_parameters->d_hierarchy = hierarchy;
   application_parameters->d_db = input_db->getDatabase("pixie3d");
   pixie3dApplication* application  = new pixie3dApplication( application_parameters );
   
   // Initialize x0
   application->setInitialConditions(0.0);

   // Create time integrator
   algs::TimeIntegratorParameters* integration_parameters = new algs::TimeIntegratorParameters();
   tbox::Pointer< tbox::Database > integrator_db = input_db->getDatabase("PredictorCorrector");
   integration_parameters->d_db = integrator_db;
   integration_parameters->dt_method = 1;
   integration_parameters->d_ic_vector = application->get_x();
   integration_parameters->d_application_strategy = application;
   integration_parameters->d_integrator_name = "PredectorCorrector";
   algs::ForwardEulerTimeIntegrator* integrator = new algs::ForwardEulerTimeIntegrator(integration_parameters);

   // Create x(t)
   tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > x_t;
   x_t = integrator->getCurrentSolution();

   // Initialize the writer
   appu::VisItDataWriter<NDIM>* d_visit_writer;
   d_visit_writer = new appu::VisItDataWriter<NDIM>("test",write_path);
   string var_name;
   stringstream stream;
   for (int i=0; i<NVAR; i++) {
      stream << "x(" << i << ")"; 
      var_name = stream.str();
      stream.str("");
      switch (i) {
         case 0:
            var_name += " - Rho";  break;
         case 1:
            var_name += " - P^1";  break;
         case 2:
            var_name += " - P^2";  break;
         case 3:
            var_name += " - P^3";  break;
         case 4:
            var_name += " - B^1";  break;
         case 5:
            var_name += " - B^2";  break;
         case 6:
            var_name += " - B^3";  break;
         case 7:
            var_name += " - Temp"; break;
      }
      const int x_id = x_t->getComponentDescriptorIndex(i);
      d_visit_writer->registerPlotQuantity(var_name,"SCALAR",x_id);
   }
   

   // Write the data
   d_visit_writer->writePlotData(hierarchy,0,0);

   // Loop through time
   double dt = integrator->getCurrentDt();
   double time = integrator->getCurrentTime();
   tbox::pout << "dt = " << dt << "\n";
   int i = 1;
   int retval;
   double last_save = 0.0;
   char buffer [50];
   while ( time < integrator->getFinalTime() ) {
      // Advance the solution
      if (i==1)
         integrator->advanceSolution(dt,true);
      else
         integrator->advanceSolution(dt,false);
      // Get the current solution and time
      x_t = integrator->getCurrentSolution();
      time = integrator->getCurrentTime();
      // Update the timestep
      dt = integrator->getNextDt(true,retval);
      // Write data
      if ( i%5==0 ) {
         sprintf(buffer,"t = %8.4f, dt = %8.6f\n",time,dt);
         tbox::pout << buffer;
      }
      if ( time >= last_save+dt_save ) {
         d_visit_writer->writePlotData(hierarchy,i,time);
         last_save += dt_save;
      }
      i++;
   }

   // That's all, folks!
   delete application;
   PetscFinalize();
   tbox::SAMRAIManager::shutdown();
   tbox::MPI::finalize();

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

