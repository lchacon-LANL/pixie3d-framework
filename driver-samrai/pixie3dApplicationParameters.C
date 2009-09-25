//
// $Id: pixie3dApplicationParameters.C 2232 2006-01-10 20:28:03Z pernice $
// $Revision: 2232 $
// $Date: 2006-01-10 13:28:03 -0700 (Tue, 10 Jan 2006) $
//
// File:  pixie3dApplicationParameters.C
// Copyright:  (c) 2005 The Regents of the University of California
// Description:  Encapasulation of parameters needed to initialize an 
//               application.
//

#include "pixie3dApplicationParameters.h"

pixie3dApplicationParameters::pixie3dApplicationParameters()
{
   d_hierarchy.setNull();
}

pixie3dApplicationParameters::pixie3dApplicationParameters(const tbox::Pointer<tbox::Database> &database) 
{
   d_hierarchy.setNull();
}

pixie3dApplicationParameters::~pixie3dApplicationParameters()
{
}


