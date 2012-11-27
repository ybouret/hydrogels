#ifndef LIBRARY_INCLUDED
#define LIBRARY_INCLUDED 1

#include "yocto/aqueous/lua.hpp"
#include "yocto/exception.hpp"

using namespace yocto;
using namespace aqueous;

//! extraneous data for each species
class SpeciesData
{
public:
    SpeciesData( double DiffusionCoeff ) throw();
    ~SpeciesData() throw();
    
    const double D;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(SpeciesData);
};


//! auto-initalized library
class Library : public library
{
public:
    explicit Library( lua_State *L );
    virtual ~Library() throw();
    
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Library);
    
};

#endif
