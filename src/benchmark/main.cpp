#include "yocto/string/vfs-utils.hpp"
#include "yocto/exception.hpp"
#include "yocto/spade/workspace.hpp"
#include "yocto/spade/rmesh.hpp"
#include "yocto/math/types.hpp"
#include "yocto/code/utils.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/fcn/functions.hpp"

#include <iostream>

using namespace yocto;
using namespace spade;
using namespace math;

typedef double        Real;
typedef array1D<Real> Array;
typedef layout1D      Layout;
typedef workspace<Layout,rmesh,Real> Workspace;


static const Real Kw = 1e-14;
const double      h0 = 1e-2;
const double      h1 = 1e-10;



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
        break;
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
    
    
    
    Real Kernel( Real t, Real x ) const
    {
        const Real args = (x/Sqrt(4*Dh*t));
        const Real d0   = h0 - Kw/h0;
        const Real d1   = h1 - Kw/h1;
        const Real d    = d1 + (d0-d1) * qerfc(args);
        return 0.5* ( d + Sqrt( d*d + 4.0 * Kw ) );
    }
    
    
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
    
    
    
    
    //--------------------------------------------------------------------------
    //-- Compute fluxes, valid for any scheme
    //--------------------------------------------------------------------------
    void compute_fluxes() throw()
    {
        for(unit_t i=0;i<imax;++i)
        {
            const unit_t ip = i+1;
            Fh[i] = -Dh * (h[ip] - h[i]) / dx;
            Fw[i] = -Dw * (w[ip] - w[i]) / dx;
        }
    }
    
    //--------------------------------------------------------------------------
    //-- Compute explicit fluxes
    //--------------------------------------------------------------------------
    void compute_explicit_increases( const double dt ) throw()
    {
        for(unit_t i=1;i<imax;++i)
        {
            const unit_t im = i-1;
            Ih[i] = - dt * (Fh[i] - Fh[im])/dx;
            Iw[i] = - dt * (Fw[i] - Fw[im])/dx;
        }
        
    }
    
    //--------------------------------------------------------------------------
    //-- Compute constrained fluxes
    //--------------------------------------------------------------------------
    void compute_constrained_increases( const double dt ) throw()
    {
        for(unit_t i=1;i<imax;++i)
        {
            const unit_t im = i-1;
            const Real   _h = h[i];
            const Real   _w = w[i];
            const Real   _s = _h+_w;
            const Real   gh = - dt * (Fh[i] - Fh[im])/dx;
            const Real   gw = - dt * (Fw[i] - Fw[im])/dx;
            const Real   d  = (gh - gw)/_s;
            
            Ih[i] =   _h * d;
            Iw[i] = - _w * d;
        }
    }
    
    
    
    
    //--------------------------------------------------------------------------
    //-- Updated explicit step with scaling
    //--------------------------------------------------------------------------
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
    
    //--------------------------------------------------------------------------
    //-- Updated constraint step with scaling
    //--------------------------------------------------------------------------
    void update_constrained( const Real scaling ) throw()
    {
        update_explicit(scaling);
        for(unit_t i=imax;i>0;--i)
            normalize(h[i], w[i], ftol);
    }
    
    
    
    //--------------------------------------------------------------------------
    //-- initialize to a pH step
    //--------------------------------------------------------------------------
    void initialize() throw()
    {
        h[0] = h0;
        w[0] = Kw/h[0];
        
        h[1] = h1;
        w[1] = Kw/h[1];
        for( unit_t i=2;i<=imax;++i)
        {
            h[i] = h[1];
            w[i] = w[1];
        }
    }
    
    void find( Real &scaling, const Array &u, const Array &du ) throw()
    {
        for( unit_t i=1; i < imax; ++i )
        {
            const Real dU = du[i];
            const Real U  = u[i];
            if( dU < -U )
            {
                //std::cerr << "u[" << i << "]=" << U << " / dU=" << dU << std::endl;
                //abort();
                const Real factor = -U/dU; assert(factor>0);
                if( scaling < 0 )
                    scaling = factor;
                else
                {
                    if(factor<scaling)
                        scaling = factor;
                }
            }
        }
    }
    
    void step(Real t,
              Real dt,
              void (Simulation::*compute_increases)(Real),
              void (Simulation::*update_fields)(Real) )
    {
        assert(compute_increases);
        assert(update_fields);
        
        Real t_end = t + dt;
        while(t<t_end)
        {
            dt = t_end - t;
            compute_fluxes();
            ((*this).*(compute_increases))(dt);
            Real scaling = -1;
            find(scaling,h,Ih);
            find(scaling,w,Iw);
            if( scaling > 0 )
            {
                std::cerr << "scaling=" << scaling << std::endl;
                scaling *= 0.5;
                dt *= scaling;
                t += dt;
                ((*this).*(update_fields))(scaling);
                continue;
            }
            else
            {
                std::cerr << "full" << std::endl;
                ((*this).*(update_fields))(1);
                break; // was a full step
            }
        }
    }
    
    void step_explicit(Real t, Real dt)
    {
        step(t,dt,&Simulation::compute_explicit_increases,&Simulation::update_explicit);
    }
    
    void step_constrained(Real t, Real dt )
    {
        step(t,dt,&Simulation::compute_constrained_increases,&Simulation::update_constrained);
    }
    
    
    void save_profile( const string &fn, Real t ) const
    {
        ios::ocstream fp( fn, false);
        for(unit_t i=0; i <= imax; ++i )
        {
            fp("%g %g %g %g\n", X[i], h[i], -log10(h[i]), -log10( Kernel(t,X[i])));
        }
    }
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Simulation);
};


static inline
double dt_round( double dt_max )
{
    //std::cerr << "dt_max=" << dt_max << std::endl;
    const double dt_log = floor( log10(dt_max) );
    //std::cerr << "dt_log=" << dt_log << std::endl;
    const double dt_one = floor(dt_max * pow(10.0,-dt_log));
    //std::cerr << "dt_one=" << dt_one << std::endl;
    return dt_one * pow(10.0,dt_log);
}


int main(int argc, char *argv[])
{
    const char *progname = _vfs::get_base_name( argv[0]);
    try
    {
        const Real alpha = 0.1;
        Simulation sim(1000,1e-2);
        const Real dt_max = alpha * (sim.dx*sim.dx) / max_of(sim.Dh,sim.Dw);
        const Real dt     = dt_round(dt_max);
        std::cerr << "dt=" << dt << std::endl;
        const Real t_run  = 4;
        Real       t      = 0;
        size_t iter_max = 1+ceil(t_run/dt);
        std::cerr << "iter_max=" << iter_max << std::endl;
        
        sim.initialize();
        sim.save_profile("h0.dat",0);
        for(size_t iter=1;iter<=iter_max;++iter)
        {
            std::cerr << "iter=" << iter << std::endl;
            sim.step_constrained(t,dt);
            t = iter*dt;
        }
        sim.save_profile("h1.dat",t);
        
        
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
