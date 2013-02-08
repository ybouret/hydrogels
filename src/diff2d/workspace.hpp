#ifndef WORKSPACE_INCLUDED
#define WORKSPACE_INCLUDED 1

#include "param.hpp"
#include "yocto/spade/variables.hpp"

typedef array1D<double> Axis;

class Workspace : public WorkspaceBase
{
public:
    explicit Workspace( const Parameters &param , library &l);
    virtual ~Workspace() throw();
    const Axis &X;
    const Axis &Y;
    Array      &porosity;
    
    void loadC( array<double> &C, const Coord &u ) const;
    void saveC( const array<double> &C, const Coord &u );
    
    variables      var;
    variables      dvar;
    
private:
    linear_handles handles;
    linear_handles handles_dC;
    YOCTO_DISABLE_COPY_AND_ASSIGN(Workspace);
};


#endif
