
#include "PCDensityRefinePatchStrategy.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/PatchLevel.h"

#include "source/AMRUtilities.h"

namespace SAMRAI {

namespace SAMRSolvers {

PCDensityRefinePatchStrategy::PCDensityRefinePatchStrategy(
      const tbox::Dimension & dim,
      BoundaryConditionParameters *parameters)
   :xfer::RefinePatchStrategy(dim)
{
   d_data_id = -1;
   d_data_idx = -1;
   d_extrapolation_order = 1;
   d_debug_print_info_level   = 0;
   d_bdry_types = new int[2*dim.getValue()];

   d_hierarchy = parameters->d_hierarchy;

   getFromInput(parameters->d_db);
}
 
PCDensityRefinePatchStrategy::~PCDensityRefinePatchStrategy()
{
  delete [] d_bdry_types;
  d_bdry_types = NULL;
}

void
PCDensityRefinePatchStrategy::getFromInput(const tbox::Pointer<tbox::Database> &db)
{
   if (db->keyExists("extrapolation_order")) {
      d_extrapolation_order = db->getInteger("extrapolation_order");
   } else {
      TBOX_ERROR("PCDensityRefinePatchStrategy"
            << " -- Required key `extrapolation_order'"
            << " missing in input.");
   }

   if (db->keyExists("boundary_conditions")) 
   {
      // get the database object for boundary conditions
      db->getIntegerArray("boundary_conditions", d_bdry_types, 2*getDim().getValue());

      for (int d=0; d<getDim().getValue(); d++) 
      {
         //	 tbox::pout << "Boundary conditions in direction " << d << " are " << d_bdry_types[2*d] << ", " << d_bdry_types[2*d+1] << std::endl; 	  
      }
   } 
   else 
   {
      TBOX_ERROR("PCDensityRefinePatchStrategy" 
            << " -- Required key `boundary_conditions'"
            << " missing in input.");
   }
}

void 
PCDensityRefinePatchStrategy::setBoundaryTypes(int* bdry_types)
{
   const int dim = getDim().getValue();

   TBOX_ASSERT(d_bdry_types!=NULL);

   for (int k=0; k<2*dim; k++) 
   {
      d_bdry_types[k] = bdry_types[k];
   }
}

void PCDensityRefinePatchStrategy::setExtrapolationOrder(int order)
{
   d_extrapolation_order = order;
}

void PCDensityRefinePatchStrategy::setPhysicalBoundaryConditions(
                                  hier::Patch& patch,
                                  const double time,
                                  const hier::IntVector& ghost_width_to_fill)
{
   if(d_extrapolation_order>0)
   {
      if(patch.inHierarchy())
      {
         if(d_debug_print_info_level>4)
         {
            tbox::pout << "Patch IS in hierarchy, level = " << patch.getPatchLevelNumber() << ", patch number = " << patch.getGlobalId() << std::endl;
         }
      }
      else
      {
         if(d_debug_print_info_level>4)
         {
            tbox::pout << "Patch NOT in hierarchy" << std::endl;
         }
      }

      /* 
       * Grab data to operate on. 
       */
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(patch.checkAllocated(d_data_id));
#endif 
      tbox::Pointer< pdat::CellData<double> > u;

      if(patch.checkAllocated(d_data_id))
      {
         u = patch.getPatchData(d_data_id);
         if(d_debug_print_info_level>4)
         {
            tbox:: pout << d_data_id << " is allocated " << std::endl;
         }
      }
      else
      {
         if(d_debug_print_info_level>4)
         {
            tbox:: pout << "ERROR:: " << d_data_id << " is NOT allocated " << std::endl;
         }
         abort();
      }

      tbox::Pointer< pdat::FaceData<double> > d;

      hier::IntVector gcw = u->getGhostCellWidth();

      const hier::Index ifirst = patch.getBox().lower();
      const hier::Index ilast =  patch.getBox().upper();

      /*
       * Determine boxes that touch the physical boundary.
       */

      const tbox::Pointer<geom::CartesianPatchGeometry > patch_geom = 
         patch.getPatchGeometry();
      const double* dx = patch_geom->getDx();
      const tbox::Array<hier::BoundaryBox > boundary = patch_geom->getCodimensionBoundaries(1);

      /* 
       * Walk the list of boxes that describe the boundary, and apply
       * boundary conditions to each segment.  Make sure we get the 
       * corners!  
       */

      int face;
      int type;

      const int dim = getDim().getValue();

      int *ilo = new int[dim];
      int *ihi = new int[dim];

      for (int axis=0; axis<dim; axis++) 
      {
         for(int side=0; side<2; side++)
         {
            hier::BoxContainer physical;
            AMRUtilities::getPhysicalBoxes(physical, patch, axis, side);
            for (hier::BoxContainerIterator l=physical.begin(); l!=physical.end(); l++) 
            {
               hier::Index pfirst = l().lower();
               hier::Index plast = l().upper();

               for(int i=0; i<dim; i++)
               {
                  ilo[i] = pfirst[i];
                  ihi[i] = plast[i];
               }

               face = 2*axis+side;
               type = 1;

               if(dim==1)
               {
               }
               else if(dim==2)
               {
               }
               else if(dim==3)
               {
               }
            }
         }
      }

      /* 
       * Next get the damned corners!
       */
      if (dim==2)
      {
         const tbox::Array<hier::BoundaryBox > corners = patch_geom->getNodeBoundaries();

         for (int i=0; i<corners.getSize(); i++) 
         {
            hier::Box bbox = corners[i].getBox();
            hier::Index pfirst = bbox.lower();
            hier::Index plast = bbox.upper();

            for(int j=0; j<dim; j++)
            {
               ilo[j] = pfirst[j];
               ihi[j] = plast[j];
            }

            face = corners[i].getLocationIndex();
            type = corners[i].getBoundaryType();

            if(d_debug_print_info_level>4)
            {
               tbox::pout << "Calling corner sets" << std::endl;
               tbox::pout << "Patch Box::" << patch.getBox() << std::endl;
               tbox::pout << "PatchEdge Box::" << bbox << std::endl;
               tbox::pout << "Face :: " << face << std::endl;
               tbox::pout << "Type :: " << type << std::endl;
               tbox::pout << "Extrapolation order " << d_extrapolation_order << std::endl;
            }

         }
      }
      else if (dim==3)
      {
         const tbox::Array<hier::BoundaryBox > edges = patch_geom->getEdgeBoundaries();

         for (int i=0; i<edges.getSize(); i++) 
         {
            hier::Box bbox = edges[i].getBox();
            hier::Index pfirst = bbox.lower();
            hier::Index plast = bbox.upper();

            for(int j=0; j<dim; j++)
            {
               ilo[j] = pfirst[j];
               ihi[j] = plast[j];
            }
            face = edges[i].getLocationIndex();
            type = edges[i].getBoundaryType();

         }

         const tbox::Array<hier::BoundaryBox > corners = patch_geom->getNodeBoundaries();

         for (int i=0; i<corners.getSize(); i++) 
         {
            hier::Box bbox = corners[i].getBox();
            hier::Index pfirst = bbox.lower();
            hier::Index plast = bbox.upper();

            for(int j=0; j<dim; j++)
            {
               ilo[j] = pfirst[j];
               ihi[j] = plast[j];
            }
            face = corners[i].getLocationIndex();
            type = corners[i].getBoundaryType();

         }
      }

      delete [] ilo;
      delete [] ihi;
   }
}

/*
************************************************************************
*                                                                      *
* Extrapolate ghost cells that are patch corners at physical           *
* boundaries.                                                          *
*                                                                      *
************************************************************************
*/
void PCDensityRefinePatchStrategy::extrapolateCornerGhostCells(
        tbox::Pointer<hier::PatchHierarchy > hierarchy,
        const int ln,
        const int var_id,
        const int var_idx)
{
   const int dim = getDim().getValue();
   int *ilo = new int[dim];
   int *ihi = new int[dim];

   tbox::Pointer<hier::PatchLevel > level = hierarchy->getPatchLevel(ln);

   for (hier::PatchLevel::Iterator p(level); p; p++) 
   {
      tbox::Pointer<hier::Patch > patch = *p;  

      const hier::Index ifirst = patch->getBox().lower();
      const hier::Index ilast =  patch->getBox().upper();
      tbox::Pointer< pdat::CellData<double> > data = 
         patch->getPatchData(var_id);

      hier::IntVector gcw = data->getGhostCellWidth();

      const tbox::Pointer<geom::CartesianPatchGeometry > patch_geometry = patch->getPatchGeometry();

      if(dim==2)
      {
         const tbox::Array<hier::BoundaryBox > corners = patch_geometry->getNodeBoundaries();
         const double* dx = patch_geometry->getDx();

         for (int i=0; i<corners.getSize(); i++)
         {
            hier::Box bbox = patch_geometry->getBoundaryFillBox(corners[i],
                  patch->getBox(),
                  gcw);
            hier::Index pfirst = bbox.lower();
            hier::Index plast = bbox.upper();

            for(int j=0; j<dim; j++)
            {
               ilo[j] = pfirst[j];
               ihi[j] = plast[j];
            }
            int face = corners[i].getLocationIndex();
            int type = corners[i].getBoundaryType();

         }      
      }
      else if(dim==3)
      {
         const tbox::Array<hier::BoundaryBox > edges =  patch_geometry->getEdgeBoundaries();
         const double* dx = patch_geometry->getDx();

         for (int i=0; i<edges.getSize(); i++) 
	   {
	     hier::Box bbox = patch_geometry->getBoundaryFillBox(edges[i],
								 patch->getBox(),
								 gcw);
	     hier::Index pfirst = bbox.lower();
	     hier::Index plast = bbox.upper();
	     
	     for(int j=0; j<dim; j++)
	       {
		 ilo[j] = pfirst[j];
		 ihi[j] = plast[j];
	       }
	     
	     int face = edges[i].getLocationIndex();
	     int type = edges[i].getBoundaryType();
	     
	   }      
	 
         const tbox::Array<hier::BoundaryBox > corners =  patch_geometry->getNodeBoundaries();

         for (int i=0; i<corners.getSize(); i++)
         {
            hier::Box bbox = patch_geometry->getBoundaryFillBox(corners[i],
                  patch->getBox(),
                  gcw);
            hier::Index pfirst = bbox.lower();
            hier::Index plast = bbox.upper();

            for(int j=0; j<dim; j++)
            {
               ilo[j] = pfirst[j];
               ihi[j] = plast[j];
            }
            int face = corners[i].getLocationIndex();
            int type = corners[i].getBoundaryType();

         }      
      }
   }

   delete [] ilo;
   delete [] ihi;

}

void 
PCDensityRefinePatchStrategy::postprocessRefine(hier::Patch& fine,
						const hier::Patch& coarse,
						const hier::Box& fine_box,
						const hier::IntVector& ratio) 
{
}

bool PCDensityRefinePatchStrategy::contains(
					    const SAMRAI::hier::Box & bigBox,
					    const SAMRAI::hier::Box & smallBox,
					    const SAMRAI::hier::IntVector & ghosts)
{
   const int dim = bigBox.getDim().getValue();

   for(int d = 0; d < dim; d++)
   {
      SAMRAI::hier::Box box = bigBox;
      if(box.lower(d) != box.upper(d))
      {
         box.grow(d,ghosts(d));
      }
      if(box.contains(smallBox))
      {
         return true;
      }
   }

   return false;
}

}
}
