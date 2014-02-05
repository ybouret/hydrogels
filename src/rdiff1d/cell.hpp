#ifndef CELL_INCLUDED
#define CELL_INCLUDED 1

#include "params.hpp"
#include "initializer.hpp"
#include "yocto/spade/rmesh.hpp"
#include "yocto/math/dat/interpolate.hpp"
using namespace math;

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
    Initializer  ini_left;
    Initializer  ini_core;
    const bool   right_wall;
    Initializer  ini_right;
    
    solution on_left;
    solution in_core;
    solution on_right;
    
    const Array1D &X;
    const Array1D &idX;
    Array1D       &pH;
    Array1D       &Q;

    vector<Array1D *> pConc;
    vector<Array1D *> pFlux;
    vector<Array1D *> pIncr;
    vector<double>    D;
    vector<int>       Z;
    
    const bool        search_front;
    const double      search_value;
    const string      search_field;
    
    // fill side/core
    void init_all(void) throw();
    
    // normalize all volumes, compute pH
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
    
    bool find_front(double &pos ) const;
    bool find_inflection(double &xx, double &yy) const;
    
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Cell);
    double weight1;
    double weight2;
    mutable interpolator<double>      front_ip;
    mutable numeric<double>::function front_fn;
    double zfront(double) throw();
};


#endif
