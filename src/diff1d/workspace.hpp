#ifndef WORKSPACE_INCLUDED
#define WORKSPACE_INCLUDED 1

#include "param.hpp"

typedef array1D<double> Axis;

class Workspace : public WorkspaceBase
{
public:
    explicit Workspace( const Parameters &param , library &l);
    virtual ~Workspace() throw();
    const Axis &X;
    
    void loadC( array<double> &C, unit_t x ) const;
    void saveC( const array<double> &C, unit_t x );
    
private:
    linear_handles handles;
    YOCTO_DISABLE_COPY_AND_ASSIGN(Workspace);
};


#endif
