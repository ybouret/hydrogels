#ifndef CELL_INCLUDED
#define CELL_INCLUDED 1


#include "initializer.hpp"
#include "workspace.hpp"
#include "yocto/spade/vtk/writer.hpp"



class Cell : public Library, public ChemSys, public Parameters, public Workspace
{
public:
    
    explicit Cell( lua_State *L);
    virtual ~Cell() throw();
    
    Initializer ini_bath;
    Initializer ini_skin;
    
    vtk_writer  vtk;
    
    void compute_fluxes();
    void compute_increases(double t, double dt);
    void update_fields(double t);
    
    double shrink;

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Cell);

};

#endif
