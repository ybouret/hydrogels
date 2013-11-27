#ifndef INIT_INCLUDED
#define INIT_INCLUDED 1

#include "chemsys.hpp"
#include "yocto/chemical/initializer.hpp"


class Initializer : public initializer
{
public:
    explicit Initializer(lua_State *L, const char *name, const collection &lib);
    virtual ~Initializer() throw();
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Initializer);
};

#endif

