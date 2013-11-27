#include "chemsys.hpp"

ChemSys:: ~ChemSys() throw()
{
    
}

ChemSys :: ChemSys( lua_State *L, const collection &lib ) : equilibria()
{
    _lua::load(L, lib, *this, "eqs" );
}


