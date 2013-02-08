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
volumes( __check_volumes(Lua::Config::Get<LUA_NUMBER>(L, "xvolumes")), __check_volumes(Lua::Config::Get<LUA_NUMBER>(L, "yvolumes")) ),
vlayout( Coord(0,0), volumes),
vtop(vlayout.upper.x-1,vlayout.upper.y-1),
length( __check_length( Lua::Config::Get<LUA_NUMBER>(L,"Lx") ), __check_length( Lua::Config::Get<LUA_NUMBER>(L,"Ly")) ),
delta( length.x/vlayout.upper.x, length.y/vlayout.upper.y),
twodel( 2.0*delta ),
noghost(),
fields( 5+l.size() * 3 )
{
    std::cerr << "*** Registering fields" << std::endl;
    //--------------------------------------------------------------------------
    // Register Fields
    //--------------------------------------------------------------------------
    for( library::const_iterator i = l.begin(); i != l.end(); ++i )
    {
        const species &sp  = **i;
        const string   spf = sp.name + "F";
        const string   spi = sp.name + "I";
        Y_SPADE_FIELD(fields, sp.name.c_str(),  Array);
        Y_SPADE_FIELD(fields, spf.c_str(),      VertexArray);
        Y_SPADE_FIELD(fields, spi.c_str(),      Array);
    }
    
    Y_SPADE_FIELD(fields, "porosity", Array);
}



Parameters:: ~Parameters() throw()
{
}
