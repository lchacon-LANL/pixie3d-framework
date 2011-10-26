//
// $Id: test_gradient.h 1832 2005-08-23 21:10:00Z bphilip $
// $Revision: 1832 $
// $Date: 2005-08-23 15:10:00 -0600 (Tue, 23 Aug 2005) $
//

/**
 * \file test_gradient.h
 * 
 * This simple example illustrates the computation of a gradient on an
 * AMR hierarchy.  
 */

#include <string>

#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"

/**
 * Extract name of input and log files from command line.
 * \param argc          Number of command line arguments.
 * \param argv          Pointers to command line arguments.
 * \param input_file    Name of input file as it appears on command line.
 * \param log_file      Name of log file as it appears on command line.
 */
void processCommandLine( int argc, char* argv[], std::string& input_file, std::string& log_file);


