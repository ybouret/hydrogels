#ifndef CHEMSYS_INCLUDED
#define CHEMSYS_INCLUDED 1

#include "library.hpp"
#include "yocto/aqueous/chemsys.hpp"

class ChemSys : public chemsys
{
public:
    explicit ChemSys( const library &l, lua_State *L );
    virtual ~ChemSys() throw();
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(ChemSys);
};

#endif
