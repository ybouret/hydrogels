#include "worker.hpp"

Worker:: ~Worker() throw() {}

Worker:: Worker(const library &l,
                lua_State     *L,
                const unit_t   s,
                const size_t   n,
                const unit_t   i_shift,
                const size_t   i_count
                ) :
ChemSys(l,L),
valid(true),
shrink(-1),
flux_shift(s),
flux_count(n),
incr_shift(i_shift),
incr_count(i_count)
{
    std::cerr << "Worker: computing flux from " << flux_shift << " to " << flux_shift+flux_count-1 << std::endl;
    std::cerr << "        computing incr from " << incr_shift << " to " << incr_shift+incr_count-1 << std::endl;
}

void Worker:: compute_fluxes( const double dx ) const throw()
{
    const Library::const_iterator end = lib.end();
    unit_t i = flux_shift;
    for(size_t n=flux_count;n>0;--n,++i)
    {
        for( library::const_iterator j=lib.begin();j != end;++j)
        {
            const species &sp = **j;
            SpeciesData   &sd = sp.get<SpeciesData>();
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

void Worker:: compute_increases(const double dt, const double dx) const throw()
{
    const Library::const_iterator end = lib.end();
    unit_t i = incr_shift;
    for(size_t n=incr_count;n>0;--n,++i)
    {
        for( library::const_iterator j=lib.begin();j != end;++j)
        {
            const species &sp = **j;
            SpeciesData   &sd = sp.get<SpeciesData>();
            assert(sd.F);
            assert(sd.I);
            const Array &F = *sd.F;
            Array       &I = *sd.I;
            I[i] = -dt * (F[i] - F[i-1]) /dx;
            //if( sp.name == "H+" ) std::cerr << "d[H+][ " << i << "]=" << I[i] << std::endl;
        }
    }
}

void Worker:: find_shrink() throw()
{
    shrink = -1;
    const Library::const_iterator end = lib.end();
    unit_t i = incr_shift;
    for(size_t n=incr_count;n>0;--n,++i)
    {
        for( library::const_iterator j=lib.begin();j != end;++j)
        {
            const species &sp = **j;
            SpeciesData   &sd = sp.get<SpeciesData>();
            assert(sd.U);
            assert(sd.I);
            const Array &U = *sd.U;
            const Array &I = *sd.I;
            assert(U[i]>=0);
            
            if( I[i] < -U[i] )
            {
                // increase is too negative
                //std::cerr << "[" << sp.name << "][" << i << "]=" << U[i] << " / delta=" << I[i] << std::endl;
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
    
}


#include "cell.hpp"

void Worker:: reduce_all(Cell &cell ) throw()
{
    valid = true;
    try
    {
        unit_t i = incr_shift;
        for(size_t n=incr_count;n>0;--n,++i)
        {
            cell.loadC( C, i );
            cell.load_dC( dC, i );
            reduce(cell.t);
            cell.save_dC( dC, i );
        }
    }
    catch(...)
    {
        valid = false;
    }
}

void Worker:: update(Cell &cell) throw()
{
    const Library::const_iterator end = lib.end();
    unit_t i = incr_shift;
    for(size_t n=incr_count;n>0;--n,++i)
    {
        size_t k = 0;
        for( library::const_iterator j=lib.begin();j != end;++j)
        {
            const species &sp = **j;
            SpeciesData   &sd = sp.get<SpeciesData>();
            assert(sd.U);
            assert(sd.I);
            Array       &U = *sd.U;
            const Array &I = *sd.I;
            C[++k] = ( U[i] += I[i] );
        }
        //cell.loadC(C,i);
        normalize(cell.t);
        cell.saveC(C,i);
    }
    
}

void Worker:: partial_update( Cell &cell ) throw()
{
    const Library::const_iterator end = lib.end();
    const double fac = cell.shrink;
    unit_t i = incr_shift;
    for(size_t n=incr_count;n>0;--n,++i)
    {
        for( library::const_iterator j=lib.begin();j != end;++j)
        {
            const species &sp = **j;
            SpeciesData   &sd = sp.get<SpeciesData>();
            assert(sd.U);
            assert(sd.I);
            Array       &U = *sd.U;
            const Array &I = *sd.I;
            U[i] += fac * I[i]; assert( U[i] >= 0 );
        }
        cell.loadC(C,i);
        normalize(cell.t);
        cell.saveC(C,i);
    }
}

