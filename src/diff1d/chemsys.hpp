#ifndef CHEMSYS_INCLUDED
#define CHEMSYS_INCLUDED 1

#include "library.hpp"
#include "yocto/chemical/equilibria.hpp"

class ChemSys : public equilibria
{
public:
    explicit ChemSys( const collection &l, lua_State *L );
    virtual ~ChemSys() throw();
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(ChemSys);
};

#endif
