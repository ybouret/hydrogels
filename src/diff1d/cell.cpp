#include "cell.hpp"
#include "yocto/lua/lua-config.hpp"
#include "yocto/code/utils.hpp"

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
crew(1,0),
workers(crew.size,as_capacity),
task_reduce( this, & Cell:: ReduceCB ),
iniBulk( "ini_bulk", *this, L ),
iniCore( "ini_core", *this, L ),
chrono(),
alpha( max_of<double>(0.5, Lua::Config::Get<LUA_NUMBER>(L, "alpha")) )
{
    
    //--------------------------------------------------------------------------
    // prepare workers
    //--------------------------------------------------------------------------
    for( size_t rank=0; rank < crew.size; ++rank )
    {
        const Worker::Ptr p( new Worker(*this,L) );
        workers.push_back(p);
    }
    
    crew.dispatch< array<Worker::Ptr>, unit_t>( workers,1,volumes-1);
    for(size_t i=1; i <= crew.size; ++i )
    {
        std::cerr << "Worker #" << i << ": " << workers[i]->start << " -> " << workers[i]->final << std::endl;
    }
    
    
    
    chrono.start();
}



void Cell:: initialize()
{
    solution S(*this);
    assert( workers.size() > 0 );
    Worker  &W = *workers.front();
    ChemSys &cs = W;
    std::cerr << "-- initialize bulk" << std::endl;
    iniBulk(cs,t);
    saveC(cs.C, X.lower);
    S.get(cs.C);
    std::cerr << "bulk=" << S << std::endl;
    
    std::cerr << "-- initialize core" << std::endl;
    iniCore(cs,t);
    for(unit_t i=X.lower+1; i <= X.upper; ++i )
    {
        saveC(cs.C,i);
    }
    S.get(cs.C);
    std::cerr << "core=" << S << std::endl;
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
#if 0
    crew.run( task_reduce );
    for(size_t i=crew.size;i>0;--i)
        if( ! workers[i]->valid )
            throw exception("Invalid Composition Found!");
#else
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
#endif
}


#include "yocto/code/utils.hpp"

void Cell:: update()
{
    //crew.cycle( task_update );
    chemsys     &cs = *workers[1];
    const size_t n  = size();
    
    // update core
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
    
    // updated wall
    for( size_t j=n;j>0;--j )
    {
        SpeciesData   &sd = *specs[j];
        assert(sd.U);
        assert(sd.I);
        const Array &U = *sd.U;
        cs.C[j] = (4*U[vtop] - U[vtop-1])/3.0;
    }
    cs.normalize(t);
    for( size_t j=n;j>0;--j )
    {
        SpeciesData   &sd = *specs[j];
        Array         &U = *sd.U;
        U[volumes] = cs.C[j];
    }
    
    
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
        update();
        t += dt;
        goto STEP;
    }
    else
    {
        shrink = 1;
        update();
    }
    return chrono.query() - stamp;
}


