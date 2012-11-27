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
            I[i] = dt * (F[i] - F[i-1]) /dx;
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
        }
    }
    catch(...)
    {
        valid = false;
    }
}
