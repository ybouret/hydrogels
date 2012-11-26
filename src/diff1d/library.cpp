#include "library.hpp"
#include "yocto/lua/lua-state.hpp"

SpeciesData:: SpeciesData( double DiffusionCoeff ) throw() :
D(DiffusionCoeff)
{
    
}

SpeciesData:: ~SpeciesData() throw()
{
}



Library:: ~Library() throw()
{
}


static inline
void SpeciesDataCtor( species &sp, lua_State *L )
{
    std::cerr << "-- Loading Data for " << sp.name << std::endl;
    if( !lua_isnumber(L,-1) )
        throw exception("Missing Diffusion Coefficicent for %s", sp.name.c_str());
    sp.make<SpeciesData,double>( lua_tonumber(L, -1));
    std::cerr << "D_{" << sp.name <<  "} = " << sp.get<SpeciesData>().D << std::endl;
}

Library:: Library( lua_State *L ) : library( sizeof(SpeciesData) )
{
    _lua::species_ctor ctor( cfunctor2(SpeciesDataCtor) );
    _lua::load(L, *this, "species", &ctor);
}

