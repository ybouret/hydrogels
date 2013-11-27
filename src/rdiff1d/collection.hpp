#ifndef COLLECTION_INCLUDED
#define COLLECTION_INCLUDED 1

#include "yocto/chemical/collection.hpp"
#include "yocto/chemical/lua/io.hpp"

using namespace yocto;
using namespace chemical;

class SpeciesData
{
public:
    explicit SpeciesData( const double user_D);
    virtual ~SpeciesData() throw();
    
    const double D;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(SpeciesData);
};

class Collection : public collection
{
public:
    explicit Collection( lua_State *L );
    virtual ~Collection() throw();
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Collection);
    void OnLoad( lua_State *L, species &sp );
    
};

#endif
