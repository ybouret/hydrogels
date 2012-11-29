#include "worker.hpp"

Worker:: ~Worker() throw() {}

Worker:: Worker(Library   &l,
                lua_State *L
                ) :
ChemSys(l,L),
start(0),
count(0),
final(0),
specs(l.specs),
valid(true)
{
}



#include "cell.hpp"

void Worker:: reduce_all( Cell &cell ) throw()
{
    valid = true;
    try
    {
        const size_t n=specs.size();
        for(unit_t i=final;i>=start;--i)
        {
            for( size_t j=n;j>0;--j)
            {
                SpeciesData &sd = *specs[j];
                const Array &U = *sd.U;
                const Array &I = *sd.I;
                C[j]  = U[i];
                dC[j] = I[i];
            }
            reduce(cell.t);
            for( size_t j=n;j>0;--j)
            {
                SpeciesData &sd = *specs[j];
                Array       &U  = *sd.U;
                U[i] = C[j];
            }
        }
    }
    catch(...)
    {
        valid = false;
    }
}


