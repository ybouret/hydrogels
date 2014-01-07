#include "initializer.hpp"

Initializer:: Initializer(lua_State *L, const char *name, const collection &lib) :
boot::loader()
{
    electroneutrality(lib);
    _lua::load(L, *this, name,lib);
}

Initializer:: ~Initializer() throw()
{
}
