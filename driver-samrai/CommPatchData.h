#ifndef included_commPatchData
#define included_commPatchData

#include <vector>

// SAMRAI headers
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "boost/shared_ptr.hpp"
#include "SAMRAI/tbox/Dimension.h"

using namespace SAMRAI;

/** \class commPatchData3dApplication
 *
 * Class to store some patch data for packing and unpacking to a communication buffer
 */
class commPatchData {
    public:

        // Constructor that initializes patch data
        commPatchData( tbox::Dimension dim );

        // Constructor that initializes patch data
        commPatchData( const boost::shared_ptr<hier::Patch>& patch, const int var_id );

        // De-constructor
        ~commPatchData();

        // Copy constructor
        commPatchData(const commPatchData& rhs);

        // Copy constructor
        commPatchData& operator=(const commPatchData& rhs);

        // Function to get the patch box
        hier::Box getBox() { return d_box; }

        // Function to get the gcw
        hier::IntVector getGCW() { return d_gcw; }

        // Function to get the variable depth
        int getDepth() { return d_depth; }

        // Function to get the data
        double* getData(int i) { return d_data[i]; }

        // Function to calculate the size of an int vector required to store the data
        size_t commBufferSize();

        // Pack the data into a preallocated int array
        void putToIntBuffer(int *buffer);

        // Unpack the data from an int array
        void getFromIntBuffer(const int *buffer);

    private:
        hier::Box d_box;
        hier::IntVector d_gcw;
        int d_depth;
        std::vector<double*> d_data;
        std::vector<size_t>  d_size;
};

#endif


