#ifndef CHEMSYS_INCLUDED
#define CHEMSYS_INCLUDED 1

#include "collection.hpp"
#include "yocto/chemical/solution.hpp"

class ChemSys : public equilibria
{
public:
    virtual ~ChemSys() throw();
    explicit ChemSys( lua_State *L, const collection &lib);
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(ChemSys);
};

#endif

