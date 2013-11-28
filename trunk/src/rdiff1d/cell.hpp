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
    const size_t M; //!< Collection::size()
    Initializer  ini_side;
    Initializer  ini_core;
    
    solution on_side;
    solution in_core;
    
    const Array1D &X;
    
    vector<Array1D *> pConc;
    vector<Array1D *> pFlux;
    vector<Array1D *> pIncr;
    vector<double>    D;
    
    // fill side/core
    void init_all(void) throw();
    
    // normalize all volumes
    void norm_all(double t);
    
    // First Pass
    void compute_fluxes();
    
    // Second Pass: return shrink factor
    double compute_increases(double dt, double t);
    
    // update with factor fac
    void   update_all(double factor);
    
    // Adaptive Step
    void step(double dt, double t);
    
    void save_xy(const string &filename) const; //!< VTK xy file
    
    double dx_min() const throw();
    double D_max() const throw();
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Cell);
};


#endif
