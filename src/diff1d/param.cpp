#include "param.hpp"
#include "yocto/lua/lua-config.hpp"


static size_t __check_volumes( size_t nv )
{
    if( nv <= 0 )
        throw exception("Invalid #volumes");
    return nv;
}

static double __check_length( double length )
{
    if( length <= 0 )
        throw exception("Invalid length");
    return length;
}

Parameters:: Parameters( const library &l, lua_State *L) :
volumes( __check_volumes(Lua::Config::Get<LUA_NUMBER>(L, "volumes")) ),
vlayout(0,volumes),
vtop(vlayout.upper-1),
length( __check_length( Lua::Config::Get<LUA_NUMBER>(L,"length") ) ),
dx( length/volumes ),
noghost(),
fields( l.size() * 3 )
{
    //--------------------------------------------------------------------------
    // Register Fields
    //--------------------------------------------------------------------------
    for( library::const_iterator i = l.begin(); i != l.end(); ++i )
    {
        const species &sp = **i;
        const string   spf = sp.name + "F";
        const string   spi = sp.name + "I";
        Y_SPADE_FIELD(fields, sp.name.c_str(),  Array);
        Y_SPADE_FIELD(fields, spf.c_str(),      Array);
        Y_SPADE_FIELD(fields, spi.c_str(),      Array);
        
    }
}

Parameters:: ~Parameters() throw()
{
}
