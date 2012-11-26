#ifndef LIBRARY_INCLUDED
#define LIBRARY_INCLUDED 1

#include "yocto/aqueous/lua.hpp"

using namespace yocto;
using namespace aqueous;

class SpeciesData
{
public:
    SpeciesData( double DiffusionCoeff ) throw();
    ~SpeciesData() throw();
    
    const double D;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(SpeciesData);
};

class Library : public library
{
public:
    explicit Library();
    virtual ~Library() throw();
    
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Library);
    
};

#endif
