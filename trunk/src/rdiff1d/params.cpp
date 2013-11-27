#include "params.hpp"

SpeciesParam:: SpeciesParam( const species &sp )  :
name( sp.name ),
flux( name + "F"),
incr( name + "I")
{
    
}


SpeciesParam:: ~SpeciesParam() throw()
{
}

const string & SpeciesParam :: key() const throw() { return name; }


Params:: ~Params() throw()
{
}

#include "yocto/code/utils.hpp"
#include "yocto/lua/lua-config.hpp"

Params:: Params( lua_State *L, const collection &lib ) :
volumes( max_of<double>(3,Lua::Config::Get<lua_Number>(L,"volumes")) ),
vlayout(0,volumes),
noghost(),
length( Lua::Config::Get<lua_Number>(L,"length") ),
F(0)
{
    std::cerr <<  "Layout: " << vlayout << std::endl;
    if(length<=0)
        throw exception("Negative Length");
    
    for( collection::const_iterator i=lib.begin(); i != lib.end(); ++i )
    {
        const species &sp = **i;
        const SpeciesParam::Pointer ptr( new SpeciesParam(sp) );
        if( ! insert(ptr) )
            throw exception("unexpected multiples species...");
        Y_SPADE_FIELD(F, ptr->name.c_str(), Array1D);
        Y_SPADE_FIELD(F, ptr->flux.c_str(), Array1D);
        Y_SPADE_FIELD(F, ptr->incr.c_str(), Array1D);
    }
}