#include "library.hpp"


SpeciesData:: SpeciesData( double DiffusionCoeff ) throw() :
D(DiffusionCoeff)
{
    
}


Library:: ~Library() throw()
{
}


static inline
void SpeciesDataCtor( species &sp, lua_State *L )
{
    std::cerr << "-- Loading Data for " << sp.name << std::endl;
}

Library:: Library( lua_State *L )
{
    _lua::species_ctor ctor( cfunctor2(SpeciesDataCtor) );
    _lua::load(L, *this, "species", &ctor);
}

