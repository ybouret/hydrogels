#ifndef CELL_INCLUDED
#define CELL_INCLUDED 1


#include "workspace.hpp"
#include "initializer.hpp"
#include "worker.hpp"

typedef threading::team Team;
typedef Team::task      Task;
typedef Team::context   Context;

class Cell : public Library, public Parameters, public Workspace
{
public:
    explicit Cell( lua_State *L );
    virtual ~Cell() throw();
    
    double              t;     //!< current time
    double              dt;    //!< current time step
    threading::team     crew;
    vector<Worker::Ptr> workers;
    Task                task_compute_fluxes;
    Task                task_compute_increases;
    Task                task_reduce;
    Task                task_shrink;
    Task                task_update;
    Initializer         iniBulk;
    Initializer         iniCore;
    
    void initialize();
    
    void compute_fluxes();
    void compute_increases();
    void reduce();
    bool find_shrink(double &s);
    void update();
    
private:
    
    void ComputeFluxesCB( const Context & ) throw();
    void ComputeIncreasesCB( const Context &) throw();
    void ReduceCB( const Context & ) throw();
    void UpdateCB( const Context & ) throw();
    void ShrinkCB( const Context & ) throw();
    
    YOCTO_DISABLE_COPY_AND_ASSIGN(Cell);
};

#endif
