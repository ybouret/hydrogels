#include "initializer.hpp"

Initializer:: Initializer(lua_State *L, const char *name) :
initializer()
{
    _lua::load(L, *this, name);
}

Initializer:: ~Initializer() throw()
{
}
