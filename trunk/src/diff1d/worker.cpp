#include "worker.hpp"

Worker:: ~Worker() throw() {}

Worker:: Worker(const library &l,
                lua_State     *L) :
ChemSys(l,L)
{
    
}

void Worker:: ComputeFluxes() throw()
{
    
}