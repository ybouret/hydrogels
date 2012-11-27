#ifndef WORKER_INCLUDED
#define WORKER_INCLUDED 1

#include "chemsys.hpp"
#include "yocto/shared-ptr.hpp"
#include "yocto/threading/team.hpp"

class Cell;
class Worker : public ChemSys
{
public:
    explicit Worker(const library &l,
                    lua_State     *L,
                    const unit_t   f_shift,
                    const size_t   f_count,
                    const unit_t   i_shift,
                    const size_t   i_count
                    );
    virtual ~Worker() throw();
    bool           valid;
    double         shrink;
    const unit_t   flux_shift;
    const size_t   flux_count;
    
    const unit_t   incr_shift;
    const unit_t   incr_count;
    
    typedef shared_ptr<Worker> Ptr;
    
    void compute_fluxes(const double dx) const throw();
    void compute_increases(const double dt, const double dx) const throw();
    void reduce_all( Cell &cell ) throw();
    void find_shrink() throw();
    void update(Cell &cell) throw();
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Worker);
};


#endif
