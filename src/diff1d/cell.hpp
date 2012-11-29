#ifndef CELL_INCLUDED
#define CELL_INCLUDED 1


#include "workspace.hpp"
#include "initializer.hpp"
#include "worker.hpp"
#include "yocto/wtime.hpp"

typedef threading::crew Crew;
typedef Crew::task      Task;
typedef Crew::context   Context;

class Cell : public Library, public Parameters, public Workspace
{
public:
    explicit Cell( lua_State *L );
    virtual ~Cell() throw();
    
    double               t;      //!< current time
    double               dt;     //!< current time step
    double               shrink; //!< shrink factor
    Crew                 crew;
    vector<Worker::Ptr>  workers;
    Task                 task_reduce;
    Initializer          iniBulk;
    Initializer          iniCore;
    wtime                chrono;
    const double         alpha;  //!< time scaling
    void   initialize();
    double max_dt() const;
    
    void compute_fluxes();
    void compute_increases();
    void reduce();
    bool found_shrink();
    void update();
    
    
    double step(double t0, double dt0);
    
    double locate( const double pH ) const;
    
private:
    
    void ComputeFluxesCB( const Context & ) throw();
    void ComputeIncreasesCB( const Context &) throw();
    void ReduceCB( const Context & ) throw();
    void UpdateCB( const Context & ) throw();
    void FindShrinkCB( const Context & ) throw();
    
    YOCTO_DISABLE_COPY_AND_ASSIGN(Cell);
};

#endif
