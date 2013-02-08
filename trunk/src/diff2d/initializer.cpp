#include "initializer.hpp"

Initializer:: ~Initializer() throw()
{
}

Initializer:: Initializer( const string &ini_name, const library &l, lua_State *L ) :
initializer(l)
{
    electroneutrality();
    _lua::load(L, *this, ini_name);
}

