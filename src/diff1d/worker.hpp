#ifndef WORKER_INCLUDED
#define WORKER_INCLUDED 1

#include "chemsys.hpp"
#include "yocto/shared-ptr.hpp"
#include "yocto/threading/crew.hpp"

class Cell;
class Worker : public ChemSys
{
public:
    explicit Worker(Library       &l,
                    lua_State     *L);
    virtual ~Worker() throw();
    
    unit_t start;
    size_t count;
    unit_t final;
    Specs &specs;
    bool   valid;
    
    typedef shared_ptr<Worker> Ptr;
    
    void reduce_all(Cell &) throw();
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Worker);
};


#endif
