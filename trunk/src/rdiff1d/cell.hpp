#ifndef CELL_INCLUDED
#define CELL_INCLUDED 1

#include "params.hpp"
#include "initializer.hpp"
#include "yocto/spade/rmesh.hpp"

typedef workspace<Layout,rmesh,double> Workspace;

class Cell :
public Collection,
public ChemSys,
public Params,
public Workspace
{
public:
    explicit Cell(lua_State *L);
    virtual ~Cell() throw();
    
    Initializer ini_side;
    Initializer ini_core;
    
    solution on_side;
    solution in_core;
    
    const Array1D &X;
    
    vector<Array1D *> pConc;
    vector<Array1D *> pFlux;
    vector<Array1D *> pIncr;
    vector<double>    D;
    
    // fill side/core
    void initialize(void) throw();
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Cell);
};


#endif
