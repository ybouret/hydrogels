#include "yocto/string/vfs-utils.hpp"
#include "yocto/exception.hpp"
#include "yocto/spade/workspace.hpp"
#include "yocto/spade/rmesh.hpp"
#include "yocto/math/types.hpp"

#include <iostream>

using namespace yocto;
using namespace spade;
using namespace math;

typedef double        Real;
typedef array1D<Real> Array;
typedef layout1D      Layout;
typedef workspace<Layout,rmesh,Real> Workspace;


static const Real Kw = 1e-14;

static size_t normalize( Real &h, Real &w, const Real ftol )
{
    size_t count = 0;
    for(;;)
    {
        ++count;
        const Real Gam = Kw - h*w;
        const Real Phi = -(w+h);
        const Real xi  = -Gam/Phi;
        const Real dh  = xi;
        const Real dw  = xi;
        h += dh;
        w += dw;
        if( Fabs(dh) > Fabs(ftol*h) ) continue;
        if( Fabs(dw) > Fabs(ftol*w) ) continue;
    }
    return count;
}

class Parameters
{
public:
    const size_t                volumes;
    const unit_t                imax;
    const Real                  length;
    const Real                  dx;
    double                      Dh,Dw;
    double                      ftol;
    const unit_t                imaxm1;
    const unit_t                imaxm2;
    double                      alpha_max;
    const Layout                sim_layout;
    const ghosts_setup          no_ghost;
    fields_setup<Layout>        fields;
    
    Parameters( size_t nv, Real len) :
    volumes(nv),
    imax(volumes),
    length( len ),
    dx(length/volumes),
    Dh(1e-8),
    Dw(1e-8),
    ftol(1e-5),
    imaxm1(imax-1),
    imaxm2(imax-2),
    alpha_max(0.1),
    sim_layout(0,volumes),
    no_ghost(),
    fields(8)
    {
        
        //-- register fields to be build
        Y_SPADE_FIELD(fields, "h", Array);
        Y_SPADE_FIELD(fields, "w", Array);
        Y_SPADE_FIELD(fields, "Fh", Array);
        Y_SPADE_FIELD(fields, "Fw", Array);
        Y_SPADE_FIELD(fields, "Ih", Array);
        Y_SPADE_FIELD(fields, "Iw", Array);
        
    }
    
    virtual ~Parameters() throw() {}
        
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Parameters);
};

class Simulation : public Parameters, public Workspace
{
public:
    Array &h;
    Array &w;
    Array &Fh;
    Array &Fw;
    Array &Ih;
    Array &Iw;
    Array &X;
    explicit Simulation( size_t nv, Real len ) :
    Parameters(nv,len),
    Workspace(sim_layout,fields,no_ghost),
    h(  (*this)["h" ].as<Array>() ),
    w(  (*this)["w" ].as<Array>() ),
    Fh( (*this)["Fh"].as<Array>() ),
    Fw( (*this)["Fw"].as<Array>() ),
    Ih( (*this)["Ih"].as<Array>() ),
    Iw( (*this)["Iw"].as<Array>() ),
    X( mesh.X() )
    {
        //-- initialiazing
        for(size_t i=0; i < imax; ++i )
        {
            X[i] = dx * i;
        }
        X[imax] = length;
    }
    
    virtual ~Simulation() throw()
    {
    }
    
    
    void compute_fluxes() throw()
    {
        for(unit_t i=0;i<imax;++i)
        {
            const unit_t ip = i+1;
            Fh[i] = Dh * (h[ip] - h[i]) / dx;
            Fw[i] = Dw * (w[ip] - w[i]) / dx;
        }
    }
    
    void compute_explicit_increases( const double dt ) throw()
    {
        for(unit_t i=1;i<imax;++i)
        {
            const unit_t im = i-1;
            Ih[i] = - dt * (Fh[i] - Fh[im]);
            Iw[i] = - dt * (Fw[i] - Fw[im]);
        }
        
    }
    
    void compute_constrained_increases( const double dt ) throw()
    {
        for(unit_t i=1;i<imax;++i)
        {
            const unit_t im = i-1;
            const Real   _h = h[i];
            const Real   _w = w[i];
            const Real   _s = _h+_w;
            const Real   gh = - dt * (Fh[i] - Fh[im]);
            const Real   gw = - dt * (Fw[i] - Fw[im]);
            const Real   d  = (gh - gw)/_s;
            
            Ih[i] =   _h * d;
            Iw[i] = - _w * d;
        }
        
               
    }
    
    void update_explicit( const Real scaling ) throw()
    {
        for(unit_t i=1;i<imax;++i)
        {
            h[i] += scaling * Ih[i];
            w[i] += scaling * Iw[i];
        }
        h[imax] = (4*h[imaxm1] - h[imaxm2])/3;
        w[imax] = (4*w[imaxm1] - w[imaxm2])/3;
    }
    
    void update_constraint( const Real scaling ) throw()
    {
        update_explicit(scaling);
        for(unit_t i=imax;i>0;--i)
            normalize(h[i], w[i], ftol);
    }
    
    
    void initialize() throw()
    {
        h[0] = 1e-2;
        w[0] = Kw/h[0];
        
        h[1] = 1e-10;
        w[1] = Kw/h[1];
        for( unit_t i=2;i<=imax;++i)
        {
            h[i] = h[1];
            w[i] = w[1];
        }
    }
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Simulation);
};




int main(int argc, char *argv[])
{
    const char *progname = _vfs::get_base_name( argv[0]);
    try
    {
        Simulation sim(100,1e-2);
        
        sim.initialize();
        
        return 0;
    }
    catch( const exception &e )
    {
        std::cerr << e.what() << std::endl;
        std::cerr << e.when() << std::endl;
    }
    catch(...)
    {
        std::cerr << "Unhandled Exception in " << progname << "!" << std::endl;
    }
    
    return 1;
}
