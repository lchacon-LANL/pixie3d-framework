#include "CommPatchData.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/SideData.h"


/************************************************************************
*  Default constructor                                                  *
************************************************************************/
commPatchData::commPatchData(tbox::Dimension dim):
    d_box(hier::Box(dim)),
    d_gcw(hier::IntVector(dim))
{ 
    d_depth = 0;
    d_data.clear();
    TBOX_ASSERT(dim.getValue()<=3);
}


/************************************************************************
*  Initializing constructor                                             *
************************************************************************/
commPatchData::commPatchData( const boost::shared_ptr<hier::Patch>& patch, const int var_id ):
    d_box(hier::Box(patch->getDim())),
    d_gcw(hier::IntVector(patch->getDim()))
{ 
    d_box = patch->getBox();
    hier::Index ifirst = d_box.lower();
    hier::Index ilast  = d_box.upper();
    int dim = d_box.getDim().getValue();
    hier::PatchData* pdat = patch->getPatchData(var_id).get();
    d_gcw = pdat->getGhostCellWidth();
    size_t Ng[3] = {1,1,1};
    for (int d=0; d<dim; d++)
        Ng[d] = ilast(d)-ifirst(d)+1+2*d_gcw(d);
    if ( dynamic_cast<pdat::CellData<double>*>(pdat)!=NULL ) {
        pdat::CellData<double>* cell = dynamic_cast<pdat::CellData<double>*>(pdat);
        d_depth = cell->getDepth();
        d_size = std::vector<size_t>(1,Ng[0]*Ng[1]*Ng[2]*d_depth);
        d_data = std::vector<double*>(1,NULL);
        d_data[0] = new double[d_size[0]];
        memcpy(d_data[0],cell->getPointer(),d_size[0]*sizeof(double));
    } else if ( dynamic_cast<pdat::FaceData<double>*>(pdat)!=NULL ) {
        pdat::FaceData<double>* face = dynamic_cast<pdat::FaceData<double>*>(pdat);
        d_depth = face->getDepth();
        d_size = std::vector<size_t>(dim,0);
        d_data = std::vector<double*>(dim,NULL);
        for (int d=0; d<dim; d++) {
            if ( d==0 )
                d_size[d] = (Ng[0]+1)*Ng[1]*Ng[2]*d_depth;
            else if ( d==1 )
                d_size[d] = Ng[0]*(Ng[1]+1)*Ng[2]*d_depth;
            else if ( d==2 )
                d_size[d] = Ng[0]*Ng[1]*(Ng[2]+1)*d_depth;
            else 
                TBOX_ERROR("error");
            d_data[d] = new double[d_size[d]];
            memcpy(d_data[d],face->getPointer(d),d_size[d]*sizeof(double));
        }
    } else if ( dynamic_cast<pdat::SideData<double>*>(pdat)!=NULL ) {
        pdat::SideData<double>* side = dynamic_cast<pdat::SideData<double>*>(pdat);
        d_depth = side->getDepth();
        d_size = std::vector<size_t>(dim,0);
        d_data = std::vector<double*>(dim,NULL);
        for (int d=0; d<dim; d++) {
            if ( d==0 )
                d_size[d] = (Ng[0]+1)*Ng[1]*Ng[2]*d_depth;
            else if ( d==1 )
                d_size[d] = Ng[0]*(Ng[1]+1)*Ng[2]*d_depth;
            else if ( d==2 )
                d_size[d] = Ng[0]*Ng[1]*(Ng[2]+1)*d_depth;
            else 
                TBOX_ERROR("error");
            d_data[d] = new double[d_size[d]];
            memcpy(d_data[d],side->getPointer(d),d_size[d]*sizeof(double));
        }
    } else {
        TBOX_ERROR("Data type is not supported yet");
    }
}


/************************************************************************
*  De-constructor                                                       *
************************************************************************/
commPatchData::~commPatchData() 
{ 
    for (size_t i=0; i<d_data.size(); i++)
        delete d_data[i];
}


/************************************************************************
*  Copy constructors                                                    *
************************************************************************/
commPatchData::commPatchData( const commPatchData& rhs ):
    d_box(rhs.d_box),
    d_gcw(rhs.d_gcw)
{ 
    d_box = rhs.d_box;
    d_gcw = rhs.d_gcw;
    d_depth = rhs.d_depth;
    d_size = rhs.d_size;
    d_data = std::vector<double*>(d_size.size(),NULL);
    for (size_t i=0; i<d_data.size(); i++) {
        d_data[i] = new double[d_size[i]*sizeof(double)];
        memcpy(d_data[i],rhs.d_data[i],d_size[i]*sizeof(double));
    }
}
commPatchData& commPatchData::operator=(const commPatchData &rhs)
{
    if ( this==&rhs ) 
        return *this;
    d_box = rhs.d_box;
    d_gcw = rhs.d_gcw;
    d_depth = rhs.d_depth;
    d_size = rhs.d_size;
    d_data = std::vector<double*>(d_size.size(),NULL);
    for (size_t i=0; i<d_data.size(); i++) {
        d_data[i] = new double[d_size[i]*sizeof(double)];
        memcpy(d_data[i],rhs.d_data[i],d_size[i]*sizeof(double));
    }
    return *this;
}


/************************************************************************
*  Calculate the buffer size                                            *
************************************************************************/
size_t commPatchData::commBufferSize() 
{
    size_t size = 0;
    // Compute the buffer size needed to communicate the d_box
    int dim = d_box.getDim().getValue();
    size += 1 + dim*2;
    // Add space to store d_gcw
    size += dim;
    // Add space to store d_depth
    size += 1;
    // Add space to store d_size and d_data
    size += 1 + d_size.size();
    for (size_t i=0; i<d_size.size(); i++)
        size += d_size[i]*sizeof(double)/sizeof(int);
    return size;
}


/************************************************************************
*  Pack the data into a preallocated int array                          *
************************************************************************/
void commPatchData::putToIntBuffer(int *buffer) 
{
    // Pack d_box
    int dim = d_box.getDim().getValue();
    hier::Index ifirst = d_box.lower();
    hier::Index ilast  = d_box.upper();
    buffer[0] = dim;
    for (int i=0; i<dim; i++)
        buffer[1+i] = ifirst[i];
    for (int i=0; i<dim; i++)
        buffer[1+dim+i] = ilast[i];
    // Pack d_gcw
    for (int i=0; i<dim; i++)
        buffer[1+2*dim+i] = d_gcw[i];
    // Pack d_depth
    buffer[1+3*dim] = d_depth;
    // Pack the data
    buffer[2+3*dim] = d_size.size();
    int* ptr = &buffer[3+3*dim];
    for (size_t i=0; i<d_size.size(); i++) {
        ptr[0] = static_cast<int>(d_size[i]);
        memcpy(&ptr[1],d_data[i],d_size[i]*sizeof(double));
        ptr += 1 + d_size[i]*sizeof(double)/sizeof(int);
    }
}


/************************************************************************
*  Unpack the data from an int array                                    *
************************************************************************/
void commPatchData::getFromIntBuffer(const int *buffer)
{
    for (size_t i=0; i<d_data.size(); i++)
        delete d_data[i];
    d_data.clear();
    d_size.clear();
    // Unpack d_box
    int dim = buffer[0];
    tbox::Dimension Dim(dim);
    hier::IntVector ifirst = hier::Index(hier::IntVector(Dim,&buffer[1]));
    hier::IntVector ilast = hier::Index(hier::IntVector(Dim,&buffer[1+dim]));
    d_box = hier::Box(hier::Index(ifirst),hier::Index(ilast),hier::BlockId(0));
    // Unpack d_gcw
    d_gcw = hier::IntVector(Dim,&buffer[1+2*dim]);
    // Unpack the depth
    d_depth = buffer[1+3*dim];
    // Unpack the data
    d_size.resize(buffer[2+3*dim],0);
    d_data.resize(buffer[2+3*dim],NULL);
    const int* ptr = &buffer[3+3*dim];
    for (size_t i=0; i<d_size.size(); i++) {
        d_size[i] = ptr[0];
        d_data[i] = new double[d_size[i]];
        memcpy(d_data[i],&ptr[1],d_size[i]*sizeof(double));
        ptr += 1 + d_size[i]*sizeof(double)/sizeof(int);
    }
}




