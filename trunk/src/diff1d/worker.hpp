#ifndef WORKER_INCLUDED
#define WORKER_INCLUDED 1

#include "chemsys.hpp"
#include "yocto/shared-ptr.hpp"
#include "yocto/threading/team.hpp"

class Worker : public ChemSys
{
public:
    explicit Worker(const library &l,
                    lua_State     *L
                    );
    virtual ~Worker() throw();
    
    
    typedef shared_ptr<Worker> Ptr;
    
    void ComputeFluxes() throw();
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Worker);
};


#endif
