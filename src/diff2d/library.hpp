#ifndef LIBRARY_INCLUDED
#define LIBRARY_INCLUDED 1

#include "types.hpp"
#include "yocto/aqueous/lua.hpp"
#include "yocto/sequence/vector.hpp"

using namespace yocto;
using namespace aqueous;

//! extraneous data for each species
class SpeciesData
{
public:
    SpeciesData( double DiffusionCoeff ) throw();
    ~SpeciesData() throw();
    
    const double  D;
    Array        *U;  //!< associated concentrations
    VertexArray  *F;  //!< associated fluxes
    Array        *I;  //!< associated increases
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(SpeciesData);
};

typedef vector<SpeciesData*> Specs;


//! auto-initalized library
class Library : public library
{
public:
    explicit Library( lua_State *L );
    virtual ~Library() throw();
    Specs specs;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Library);
    
};

#endif
