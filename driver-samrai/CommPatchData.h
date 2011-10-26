#ifndef included_commPatchData
#define included_commPatchData

// SAMRAI headers
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/Pointer.h"
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
        commPatchData( const tbox::Pointer<hier::Patch>& patch, const int var_id );

        // De-constructor
        ~commPatchData();

        // Copy constructor
        commPatchData(const commPatchData& rhs);

        // Function to get the patch box
        hier::Box getBox() { return box; }

        // Function to get the gcw
        hier::IntVector getGCW() { return gcw; }

        // Function to get the variable depth
        int getDepth() { return depth; }

        // Function to get the data
        double* getData() { return data; }

        // Function to calculate the size of an int vector required to store the data
        size_t commBufferSize();

        // Pack the data into a preallocated int array
        void putToIntBuffer(int *buffer);

        // Unpack the data from an int array
        void getFromIntBuffer(int *buffer);

    private:
        hier::Box box;
        hier::IntVector gcw;
        int depth;
        double *data;
        bool allocated_data;
};

#endif


