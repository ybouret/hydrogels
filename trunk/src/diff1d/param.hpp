#ifndef PARAM_INCLUDED
#define PARAM_INCLUDED 1

#include "library.hpp"

class Parameters
{
public:
    explicit Parameters(const collection &l, lua_State *L);
    virtual ~Parameters() throw();
    
    const size_t       volumes; //!< finite volumes
    const Layout       vlayout; //!< 0..volumes
    const unit_t       vtop;    //!< vlayout.upper-1
    const double       length;  //!< space is [0:length];
    const double       dx;      //!< length/volumes
    const ghosts_setup noghost;
    Fields             fields;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Parameters);
};

#endif
