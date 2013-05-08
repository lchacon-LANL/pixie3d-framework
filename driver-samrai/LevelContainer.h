#ifndef included_levelcontainer
#define included_levelcontainer

#include <vector>
#include <map>
#include "SAMRAI/tbox/Array.h"
#include "boost/shared_ptr.hpp"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/VariableContext.h"
#include "PatchContainer.h"

using namespace SAMRAI;


class LevelContainer{
public:

    //! Constructor
    /*!
     * Level container constructor.  This creates and initializes the level container.
     * Internally this creates a patch container object for each patch. 
     * @param hierarchy     Patch hierarchy
     * @param level         Patch level
     * @param n_var         The number of primary variables
     * @param u0_id         The ids of the equilibrium primary variables
     * @param u_id          The ids of the primary variables
     * @param n_auxs        The number of scalar auxillary variables
     * @param auxs_id       The ids of the scalar auxillary variables
     * @param n_auxv        The number of vector auxillary variables
     * @param auxv_id       The ids of the vector auxillary variables
     * @param flux_id       The id of the flux vector
     * @param src_id        The id of the src vector
     */
    LevelContainer( boost::shared_ptr<hier::PatchHierarchy> hierarchy, 
        boost::shared_ptr<hier::PatchLevel> level, 
        int n_var, int *u0_id, int *u_id, int n_auxs, int *auxs_id, 
        int n_auxv, int *auxv_id, int flux_id, int src_id );

    //! De-constructor
    ~LevelContainer();
    
    //! Get the pointer to a patch pointer
    void *getPtr( boost::shared_ptr<hier::Patch> patch );

private:
    // Prevent copying of the level container
    LevelContainer() {}
    // Store a pointer to the hierarchy
    boost::shared_ptr<hier::PatchHierarchy> d_hierarchy;
    // Store a pointer to the level
    boost::shared_ptr<hier::PatchLevel> d_level;
    // Store an array to translate the local patch ids to an index
    std::map<int,int> patch_map;
    // Store the ids of the data
    int n_var, *u0_id, *u_id, n_auxs, *auxs_id, n_auxv, *auxv_id, flux_id, src_id;
    // Store the pointers to the data 
    std::vector<PatchContainer*> data;
};

#endif


