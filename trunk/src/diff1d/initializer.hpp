#ifndef INITIALIZER_INCLUDED
#define INITIALIZER_INCLUDED 1

#include "yocto/chemical/initializer.hpp"
#include "chemsys.hpp"

//! initalizer with built-in electroneutrality
class Initializer : public initializer
{
public:
    explicit Initializer( const string &ini_name, const library &l, lua_State *L);
    virtual ~Initializer() throw();
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Initializer);
};

#endif

