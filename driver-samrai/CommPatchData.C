#include "CommPatchData.h"

#include "SAMRAI/pdat/CellData.h"


/************************************************************************
*  Default constructor                                                  *
************************************************************************/
commPatchData::commPatchData(tbox::Dimension dim):
    box(hier::Box(dim)),
    gcw(hier::IntVector(dim))
{ 
    depth = 0;
    data = NULL;
    allocated_data = false;
}


/************************************************************************
*  Initializing constructor                                             *
************************************************************************/
commPatchData::commPatchData( const boost::shared_ptr<hier::Patch>& patch, const int var_id ):
    box(hier::Box(patch->getDim())),
    gcw(hier::IntVector(patch->getDim()))
{ 
    box = patch->getBox();
    boost::shared_ptr<pdat::CellData<double> > pdat = 
        boost::dynamic_pointer_cast<pdat::CellData<double> >(patch->getPatchData(var_id));
    gcw = pdat->getGhostCellWidth();
    depth = pdat->getDepth();
    data = pdat->getPointer();
    allocated_data = false;
}


/************************************************************************
*  De-constructor                                                       *
************************************************************************/
commPatchData::~commPatchData() 
{ 
    if ( allocated_data && data!=NULL )
        delete [] data;
}

/************************************************************************
*  Copy constructor                                                     *
************************************************************************/
commPatchData::commPatchData( const commPatchData& rhs ):
    box(rhs.box),
    gcw(rhs.gcw)
{ 
    box = rhs.box;
    gcw = rhs.gcw;
    depth = rhs.depth;
    allocated_data = rhs.allocated_data;
    if ( allocated_data && rhs.data!=NULL) {
        int dim = box.getDim().getValue();
        hier::Index ifirst = box.lower();
        hier::Index ilast  = box.upper();
        size_t data_size = depth;
        for (int i=0; i<dim; i++)
            data_size *= ilast(i)-ifirst(i)+1+2*gcw(i);
        data = new double[data_size];
        for (size_t i=0; i<data_size; i++)
            data[i] = rhs.data[i];
    } else {
        data = rhs.data;
    }
}


/************************************************************************
*  Calculate the buffer size                                            *
************************************************************************/
size_t commPatchData::commBufferSize() 
{
    size_t size = 0;
    // Compute the buffer size needed to communicate the box
    int dim = box.getDim().getValue();
    size += 1 + dim*2;
    // Add space to store gcw
    size += dim;
    // Add space to store depth
    size += 1;
    // Add space to store data
    if ( dim==0 || depth==0 )
        return size;
    hier::Index ifirst = box.lower();
    hier::Index ilast  = box.upper();
    size_t data_size = depth;
    for (int i=0; i<dim; i++)
        data_size *= ilast(i)-ifirst(i)+1+2*gcw(i);
    size += data_size*sizeof(double)/sizeof(int);
    return size;
}


/************************************************************************
*  Pack the data into a preallocated int array                          *
************************************************************************/
void commPatchData::putToIntBuffer(int *buffer) 
{
    // Pack the box
    int dim = box.getDim().getValue();
    hier::Index ifirst = box.lower();
    hier::Index ilast  = box.upper();
    buffer[0] = dim;
    for (int i=0; i<dim; i++)
        buffer[1+i] = ifirst[i];
    for (int i=0; i<dim; i++)
        buffer[1+dim+i] = ilast[i];
    // Pack gcw
    for (int i=0; i<dim; i++)
        buffer[1+2*dim+i] = gcw[i];
    // Pack the depth
    buffer[1+3*dim] = depth;
    // Pack the data
    if ( dim==0 || depth==0 )
        return;
    if ( data==NULL ) 
        TBOX_ERROR("Internal Error");
    double *ptr = (double*) &buffer[2+3*dim];
    size_t data_size = depth;
    for (int i=0; i<dim; i++)
        data_size *= ilast(i)-ifirst(i)+1+2*gcw(i);
    for (size_t i=0; i<data_size; i++)
        ptr[i] = data[i];
}


/************************************************************************
*  Unpack the data from an int array                                    *
************************************************************************/
void commPatchData::getFromIntBuffer(int *buffer)
{
    // Unpack the box
    int dim = buffer[0];
    tbox::Dimension Dim(dim);
    hier::IntVector ifirst = hier::Index(hier::IntVector(Dim,&buffer[1]));
    hier::IntVector ilast = hier::Index(hier::IntVector(Dim,&buffer[1+dim]));
    box = hier::Box(hier::Index(ifirst),hier::Index(ilast),hier::BlockId(0));
    // Unpack gcw
    gcw = hier::IntVector(Dim,&buffer[1+2*dim]);
    // Unpack the depth
    depth = buffer[1+3*dim];
    // Unpack the data
    if ( data!=NULL && allocated_data ) 
        delete [] data;
    data = NULL;
    if ( dim==0 || depth==0 )
        return;
    size_t data_size = depth;
    for (int i=0; i<dim; i++)
        data_size *= ilast(i)-ifirst(i)+1+2*gcw(i);
    data = new double[data_size];
    double *ptr = (double*) &buffer[2+3*dim];
    for (size_t i=0; i<data_size; i++)
        data[i] = ptr[i];
}




