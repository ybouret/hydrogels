#include "collection.hpp"
#include "yocto/lua/lua.hpp"
#include "yocto/exception.hpp"


SpeciesData:: ~SpeciesData() throw()
{
    
}

SpeciesData:: SpeciesData( const double user_D ) :
D(user_D)
{
    
}

Collection:: ~Collection() throw()
{
}


Collection:: Collection( lua_State *L ) :
collection()
{
    // load from config file
    species_ctor cb( this, & Collection::OnLoad );
    _lua::load(L, *this, "species", &cb);
    
    // check everyone has data
    for( iterator i=begin(); i != end(); ++i )
    {
        const species &sp = **i;
        if(!sp.data.is_active())
            throw exception("Missing Data for '%s'", sp.name.c_str());
    }
}

void Collection:: OnLoad( lua_State *L, species &sp )
{
    std::cerr << "\t...loading parameters for '" << sp.name << "'" << std::endl;
    if(lua_gettop(L)<=0)
        throw exception("OnLoad shouldn't be called");
    const char *tid = lua_typename(L, lua_type(L,-1));
    //std::cerr << "parameter is '" << tid << "'" << std::endl;
    if(!lua_isnumber(L, -1))
        throw exception("Expecting a diffusion coefficient!");
    
    sp.data.build<SpeciesData,double>( lua_tonumber(L, -1) );
}