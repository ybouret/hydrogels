#include "cell.hpp"
#include "yocto/lua/lua-config.hpp"

Cell:: ~Cell() throw()
{
}

Cell:: Cell( lua_State *L ) :
Library(L),
Parameters(*this,L),
Workspace(*this,*this),
crew(),
workers(crew.size,as_capacity),
task_compute_fluxes( this, & Cell::ComputeFluxesCB),
task_compute_increases( this, & Cell::ComputeIncreasesCB),
task_reduce( this, & Cell:: ReduceCB),
iniBulk( "ini_bulk", *this, L ),
iniCore( "ini_core", *this, L )
{
    
    //--------------------------------------------------------------------------
    // prepare workers
    //--------------------------------------------------------------------------
    size_t flux_total = volumes;
    unit_t flux_shift = 0;
    
    size_t incr_total = volumes-1;
    unit_t incr_shift = 1;
    
    for( size_t rank=0; rank < crew.size; ++rank )
    {
        const size_t flux_count = flux_total/(crew.size-rank);
        const size_t incr_count = incr_total/(crew.size-rank);
        Worker::Ptr p( new Worker(*this,L,flux_shift,flux_count,incr_shift,incr_count) );
        workers.push_back(p);
        flux_total -= flux_count;
        flux_shift += flux_count;
        
        incr_total -= incr_count;
        incr_shift += incr_count;
    }
}


void Cell:: initialize()
{
    
    assert( workers.size() > 0 );
    Worker  &W = *workers.front();
    ChemSys &cs = W;
    std::cerr << "-- initialize bulk" << std::endl;
    iniBulk(cs,t);
    saveC(cs.C, X.lower);
    
    std::cerr << "-- initialize core" << std::endl;
    iniCore(cs,t);
    for(unit_t i=X.lower+1; i <= X.upper; ++i )
    {
        saveC(cs.C,i);
    }
}

void Cell:: ComputeFluxesCB(const Context &ctx) throw()
{
    const Worker &worker = *workers[ ctx.indx ];
    if(false)
    {
        scoped_lock guard(ctx.access);
        std::cerr << "Worker Fluxes: " << worker.flux_shift << "," << worker.flux_count << std::endl;
    }
    worker.compute_fluxes(dx);
}

void Cell:: compute_fluxes()
{
    crew.cycle( task_compute_fluxes );
}

void Cell:: ComputeIncreasesCB( const Context &ctx ) throw()
{
    const Worker &worker = *workers[ ctx.indx ];
    if(false)
    {
        scoped_lock guard(ctx.access);
        std::cerr << "Worker Increases: " << worker.incr_shift << "," << worker.incr_count << std::endl;
    }
    worker.compute_increases(dt,dx);
}

void Cell:: compute_increases()
{
    crew.cycle( task_compute_increases );
}

void Cell:: ReduceCB(const Context &ctx ) throw()
{
    Worker &worker = *workers[ ctx.indx ];
    worker.reduce_all(*this);
}

void Cell:: reduce()
{
    crew.cycle( task_reduce );
    for(size_t i=crew.size;i>0;--i)
        if( ! workers[i]->valid )
            throw exception("Invalid Composition Found!");
}
