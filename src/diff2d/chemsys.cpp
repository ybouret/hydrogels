#include "chemsys.hpp"
#include "yocto/lua/lua-config.hpp"

ChemSys:: ~ChemSys() throw()
{
    
}

ChemSys:: ChemSys( const library &l, lua_State *L ) :
chemsys(l, Lua::Config::Get<LUA_NUMBER>(L, "ftol") )
{
    std::cerr << "*** ftol=" << ftol << std::endl;
    _lua::load(L, *this, "chemsys");
    build();
}
