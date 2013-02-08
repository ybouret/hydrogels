#ifndef PARAM_INCLUDED
#define PARAM_INCLUDED 1

#include "library.hpp"

class Parameters
{
public:
    explicit Parameters(const library &l, lua_State *L);
    virtual ~Parameters() throw();
    
    const Coord        volumes; //!< #volumes
    const Layout       vlayout; //!< 0..volumes
    const Coord        vtop;    //!< vlayout.upper-1
    const Vertex       length;  //!< space is [0:length];
    const Vertex       delta;   //!< length/volumes
    const Vertex       twodel;  //!< 2*delta
    const ghosts_setup noghost;
    Fields             fields;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Parameters);
};

#endif
