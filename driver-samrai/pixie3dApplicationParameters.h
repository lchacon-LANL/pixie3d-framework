//
// $Id: pixie3dApplicationParameters.h 2232 2006-01-10 20:28:03Z pernice $
// $Revision: 2232 $
// $Date: 2006-01-10 13:28:03 -0700 (Tue, 10 Jan 2006) $
//
// File:  pixie3dApplicationParameters.h
// Copyright:  (c) 2005 The Regents of the University of California
// Description:  Encapasulation of parameters needed to initialize an 
//               application.
//

#ifndef included_pixie3d_application_parameters
#define included_pixie3d_application_parameters

#include <string>

#include "tbox/Database.h"
#include "PatchHierarchy.h"

#include "ApplicationParameters.h"

#ifndef LACKS_NAMESPACE
using namespace SAMRAI;
#endif

/** \class pixie3dApplicationParameters
 *
 * Class pixie3dApplicationParameters provides a uniform mechanism to pass
 * initialization parameters when constructing a pixie3dApplication.
 */
class pixie3dApplicationParameters : public ApplicationParameters
{
public:
   
   // Empty constructor.
   pixie3dApplicationParameters();

   // Construct and initialize a parameter list according to input
   // data.  See Application for a list of required and optional keywords.
   pixie3dApplicationParameters( const tbox::Pointer<tbox::Database> &database );

   // Destructor.
   virtual ~pixie3dApplicationParameters();

   // Computational grid where problem is solved.
   tbox::Pointer< hier::PatchHierarchy<NDIM> > d_hierarchy;

   // Database
   tbox::Pointer< tbox::Database > d_db;

};
#endif
