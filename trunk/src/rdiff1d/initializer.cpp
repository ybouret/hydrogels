#include "initializer.hpp"

Initializer:: Initializer(lua_State *L, const char *name, const collection &lib) :
initializer()
{
    electroneutrality(lib);
    _lua::load(L, *this, name);
}

Initializer:: ~Initializer() throw()
{
}
