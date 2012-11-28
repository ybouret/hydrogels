#include "cell.hpp"
#include "yocto/lua/lua-config.hpp"

Cell:: ~Cell() throw()
{
}

Cell:: Cell( lua_State *L ) :
Library(L),
Parameters(*this,L),
Workspace(*this,*this),
t(0),
dt(0),
shrink(0),
specs( size(), as_capacity ),
crew(1,0),
workers(crew.size,as_capacity),
task_compute_fluxes( this, & Cell::ComputeFluxesCB),
task_compute_increases( this, & Cell::ComputeIncreasesCB),
task_reduce( this, & Cell:: ReduceCB ),
task_find_shrink( this, & Cell:: FindShrinkCB ),
task_update( this, & Cell:: UpdateCB ),
task_partial_update( this, &Cell::PartialUpdateCB),
iniBulk( "ini_bulk", *this, L ),
iniCore( "ini_core", *this, L ),
chrono()
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
    
    for( library::iterator i = begin(); i != end(); ++i )
        specs.push_back( & ( (**i).get<SpeciesData>() ) );
    
    chrono.start();
}


double Cell:: max_dt() const
{
    double Dmax = 0;
    for( library::const_iterator i = begin(); i != end(); ++i )
    {
        const double D = (**i).get<SpeciesData>().D;
        if( D > Dmax ) Dmax = D;
    }
    if( Dmax <= 0 )
        throw exception("Invalid Diffusion Coefficients");
    return 0.1 * (dx*dx) / Dmax;
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
    worker.compute_fluxes(dx);
}

void Cell:: compute_fluxes()
{
    //crew.cycle( task_compute_fluxes );
    const size_t n = size();
    for(unit_t i=vtop;i>=0;--i)
    {
        for( size_t j=n;j>0;--j)
        {
            SpeciesData   &sd = *specs[j];
            assert(sd.U);
            assert(sd.F);
            const Array &U = *sd.U;
            Array       &F = *sd.F;
            assert(U[i+1]>=0);
            assert(U[i]>=0);
            F[i] = - sd.D * ( U[i+1] - U[i] )/dx;
        }
    }
}

void Cell:: ComputeIncreasesCB( const Context &ctx ) throw()
{
    const Worker &worker = *workers[ ctx.indx ];
    worker.compute_increases(dt,dx);
}

void Cell:: compute_increases()
{
    //crew.cycle( task_compute_increases );
    const size_t n = size();
    for( unit_t i=vtop;i>0;--i)
    {
        for( size_t j=n;j>0;--j )
        {
            SpeciesData   &sd = *specs[j];
            assert(sd.F);
            assert(sd.I);
            const Array &F = *sd.F;
            Array       &I = *sd.I;
            I[i] = -dt * (F[i] - F[i-1]) /dx;
        }
    }
    
}

void Cell:: ReduceCB(const Context &ctx ) throw()
{
    Worker &worker = *workers[ ctx.indx ];
    worker.reduce_all(*this);
}

void Cell:: reduce()
{
    /*
     crew.cycle( task_reduce );
     for(size_t i=crew.size;i>0;--i)
     if( ! workers[i]->valid )
     throw exception("Invalid Composition Found!");
     */
    chemsys     &cs = *workers[1];
    const size_t n  = size();
    for( unit_t i=vtop;i>0;--i)
    {
        for( size_t j=n;j>0;--j )
        {
            SpeciesData &sd = *specs[j];
            assert(sd.U);
            assert(sd.I);
            const Array &U = *sd.U;
            const Array &I = *sd.I;
            cs.C[j]  = U[i];
            cs.dC[j] = I[i];
        }
        cs.reduce(t);
        for( size_t j=n;j>0;--j )
        {
            SpeciesData &sd = *specs[j];
            Array       &I  = *sd.I;
            I[i] = cs.dC[j];
        }
    }
}

void Cell:: UpdateCB( const Context &ctx ) throw()
{
    workers[ctx.indx]->update(*this);
}

void Cell:: update()
{
    //crew.cycle( task_update );
    chemsys     &cs = *workers[1];
    const size_t n  = size();
    for( unit_t i=vtop;i>0;--i)
    {
        for( size_t j=n;j>0;--j )
        {
            SpeciesData   &sd = *specs[j];
            assert(sd.U);
            assert(sd.I);
            Array       &U = *sd.U;
            const Array &I = *sd.I;
            cs.C[j] = ( U[i] += I[i] );
        }
        cs.normalize(t);
        for( size_t j=n;j>0;--j )
        {
            SpeciesData   &sd = *specs[j];
            Array         &U = *sd.U;
            U[i] = cs.C[j];
        }
    }
    
    
}


void Cell:: PartialUpdateCB( const Context &ctx ) throw()
{
    workers[ctx.indx]->partial_update(*this);
}

void Cell:: partial_update()
{
    //crew.cycle( task_partial_update );
    chemsys     &cs = *workers[1];
    const size_t n  = size();
    for( unit_t i=vtop;i>0;--i)
    {
        for( size_t j=n;j>0;--j )
        {
            SpeciesData   &sd = *specs[j];
            assert(sd.U);
            assert(sd.I);
            Array       &U = *sd.U;
            const Array &I = *sd.I;
            cs.C[j] = ( U[i] += shrink * I[i] );
        }
        cs.normalize(t);
        for( size_t j=n;j>0;--j )
        {
            SpeciesData   &sd = *specs[j];
            Array         &U = *sd.U;
            U[i] = cs.C[j];
        }
    }
    
}

void Cell:: FindShrinkCB(const Context &ctx ) throw()
{
    workers[ctx.indx]->find_shrink();
    
}

bool Cell:: found_shrink()
{
    shrink = -1;
#if 0
    crew.cycle( task_find_shrink );
    for(size_t i=crew.size;i>0;--i)
    {
        const double ws = workers[i]->shrink;
        if( ws >= 0 )
        {
            if( shrink < 0 )
                shrink = ws;
            else
                if( ws < shrink ) shrink=ws;
        }
    }
#endif
    
    const size_t n = size();
    for( unit_t i=vtop;i>0;--i)
    {
        for( size_t j=n;j>0;--j)
        {
            SpeciesData   &sd = *specs[j];
            assert(sd.U);
            assert(sd.I);
            const Array &U = *sd.U;
            const Array &I = *sd.I;
            assert(U[i]>=0);
            
            if( I[i] < -U[i] )
            {
                // increase is too negative
                const double fac = -U[i] / I[i]; assert(fac>=0);
                if( shrink < 0 )
                    shrink=fac;
                else
                {
                    if( fac < shrink )
                        shrink=fac;
                }
            }
        }
    }
    
    return shrink >= 0;
}


double Cell:: step(double t0, double dt0)
{
    
    const double stamp = chrono.query();
    {
        //-- initialize
        t  = t0;
        const double t_end = t0 + dt0;
        
    STEP:
        dt = t_end - t;
        //-- forward physics
        compute_fluxes();
        compute_increases();
        
        //-- apply chemistry
        reduce();
        if( found_shrink() )
        {
            shrink *= 0.5;
            dt     *= shrink;
            partial_update();
            t += dt;
            goto STEP;
        }
        else
            update();
    }
    return chrono.query() - stamp;
}


