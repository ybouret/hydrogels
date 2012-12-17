#include "yocto/string/vfs-utils.hpp"
#include "yocto/exception.hpp"
#include "yocto/spade/workspace.hpp"
#include "yocto/spade/rmesh.hpp"
#include "yocto/math/types.hpp"
#include "yocto/code/utils.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/fcn/functions.hpp"
#include "yocto/wtime.hpp"
#include "yocto/math/ode/stiff-drvkr.hpp"
#include "yocto/eta.hpp"

#include <iostream>

using namespace yocto;
using namespace spade;
using namespace math;

typedef double                         Real;
typedef array1D<Real>                  Array;
typedef layout1D                       Layout;
typedef workspace<Layout,rmesh,Real>   Workspace;
typedef ode::stiff_drvkr<Real>::type   Solver;
typedef ode::field<Real>::type         DiffEq;
typedef ode::field<Real>::diff         Jacobn;

static const Real Kw = 1e-14;
const        Real h0 = 1e-2;
const        Real h1 = 1e-10;

const Real k2 = 1.40e11;
const Real kd = Kw*k2;

static Real get_rate( const Real h, const Real w )
{
    return kd-k2*h*w;
}

//! Newton step
/**
 \return #step
 */
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
    ftol(1e-7),
    imaxm1(imax-1),
    imaxm2(imax-2),
    sim_layout(0,volumes),
    no_ghost(),
    fields(8)
    {
        
        //-- register fields to be build
        Y_SPADE_FIELD(fields, "h",    Array);
        Y_SPADE_FIELD(fields, "w",    Array);
        Y_SPADE_FIELD(fields, "Fh",   Array);
        Y_SPADE_FIELD(fields, "Fw",   Array);
        Y_SPADE_FIELD(fields, "Ih",   Array);
        Y_SPADE_FIELD(fields, "Iw",   Array);
        Y_SPADE_FIELD(fields, "step", Array);
        
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
    Array       &h;
    Array       &w;
    Array       &Fh;
    Array       &Fw;
    Array       &Ih;
    Array       &Iw;
    Array       &step; //!< for ODE
    vector<Real> Y;
    Solver       ODE;
    DiffEq       drvs;
    Jacobn       djac;
    wtime        chrono;
    Array       &X;
    double       t_diff;
    double       t_chem;
    double       t_proj;
    
    explicit Simulation( size_t nv, Real len ) :
    Parameters(nv,len),
    Workspace(sim_layout,fields,no_ghost),
    h(    (*this)["h" ].as<Array>() ),
    w(    (*this)["w" ].as<Array>() ),
    Fh(   (*this)["Fh"].as<Array>() ),
    Fw(   (*this)["Fw"].as<Array>() ),
    Ih(   (*this)["Ih"].as<Array>() ),
    Iw(   (*this)["Iw"].as<Array>() ),
    step( (*this)["step"].as<Array>() ),
    Y(2,0.0),
    ODE(ftol),
    drvs( this, & Simulation:: diffeq ),
    djac( this, & Simulation:: jacobn ),
    X( mesh.X() ),
    t_diff(0),
    t_chem(0),
    t_proj(0)
    {
        //-- initialiazing
        for(unit_t i=0; i < imax; ++i )
        {
            X[i] = dx * i;
        }
        X[imax] = length;
        
        //--
        ODE.start(2);
    }
    
    virtual ~Simulation() throw()
    {
    }
    
    void diffeq( array<Real> &dYdt, Real t, const array<Real> &Y ) const throw()
    {
        const Real rate = get_rate(Y[1], Y[2]);
        //std::cerr << "rate=" << rate << std::endl;
        dYdt[1] = rate; //! F[1]
        dYdt[2] = rate; //! F[2]
    }
    
    void jacobn( array<Real> &dFdt, matrix<Real> &dFdY, Real t, const array<Real> &Y ) const throw()
    {
        dFdt[1] = 0;
        dFdt[2] = 0;
        
        dFdY[1][1] = dFdY[2][1] = -k2 * Y[2];
        dFdY[1][2] = dFdY[2][2] = -k2 * Y[1];
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
    //-- Compute relaxed fluxes
    //--------------------------------------------------------------------------
    void compute_relaxed_increases( const double dt ) throw()
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
    //-- Straight explicit update
    //--------------------------------------------------------------------------
    inline void add_increases( const Real scaling ) throw()
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
    //-- Updated explicit step with scaling
    //--------------------------------------------------------------------------
    inline void update_explicit( const Real scaling, const Real dt )
    {
        for(unit_t i=1;i<imax;++i)
        {
            Y[1] = ( h[i] += scaling * Ih[i] );
            Y[2] = ( w[i] += scaling * Iw[i] );
            
            const Real rate = get_rate(Y[1], Y[2]);
            if(rate<0)
            {
                const Real max_step = 0.9*min_of( Y[1], Y[2] ) /  -rate;
                step[i] = min_of(step[i],max_step);
            }
            //ODE( drvs, djac, Y, 0, dt, step[i] );
            
            h[i] = Y[1];
            w[i] = Y[2];
        }
        h[imax] = (4*h[imaxm1] - h[imaxm2])/3;
        w[imax] = (4*w[imaxm1] - w[imaxm2])/3;
    }
    
    //--------------------------------------------------------------------------
    //-- Updated constraint step with scaling
    //--------------------------------------------------------------------------
    inline void update_relaxed( const Real scaling, const Real dt) throw()
    {
        
        for( unit_t i=1;i<imax;++i)
        {
            h[i] += scaling * Ih[i];
            w[i] += scaling * Iw[i];
        }
        h[imax] = (4*h[imaxm1] - h[imaxm2])/3;
        w[imax] = (4*w[imaxm1] - w[imaxm2])/3;
        normalize(h[imax], w[imax], ftol);
    }
    
    
    //--------------------------------------------------------------------------
    //-- initialize to a pH step
    //--------------------------------------------------------------------------
    void initialize() throw()
    {
        chrono.start();
        h[0] = h0;
        w[0] = Kw/h[0];
        
        h[1] = h1;
        w[1] = Kw/h[1];
        for( unit_t i=2;i<=imax;++i)
        {
            h[i] = h[1];
            w[i] = w[1];
        }
        
        for(unit_t i=0;i<=imax;++i) step[i] = 1e-6;
        t_diff = t_chem = t_proj = 0;
    }
    
    //--------------------------------------------------------------------------
    //-- find the minimum scaling factor in case of overshoot
    //--------------------------------------------------------------------------
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
    
    void run(Real t,
             Real dt,
             void (Simulation::*compute_increases)(Real),
             void (Simulation::*update_fields)(Real,Real)
             )
    {
        assert(compute_increases);
        assert(update_fields);
        Real t_end = t + dt;
        while(t<t_end)
        {
            dt = t_end - t;
            const double t_enter = chrono.query();
            compute_fluxes();
            ((*this).*(compute_increases))(dt);
            t_diff += t_enter - chrono.query();
            Real scaling = -1;
            find(scaling,h,Ih);
            find(scaling,w,Iw);
            
            if( scaling > 0 )
            {
                //std::cerr << "scaling=" << scaling << std::endl;
                scaling *= 0.5;
                dt *= scaling;
                ((*this).*(update_fields))(scaling,dt);
                t += dt;
                continue;
            }
            else
            {
                ((*this).*(update_fields))(1,dt);
                break; // was a full step
            }
        }
        
    }
    
    
    inline void run_explicit(Real t, Real dt)
    {
        run(t,
            dt,
            &Simulation::compute_explicit_increases,
            &Simulation::update_explicit);
    }
    
    inline void run_relaxed(Real t, Real dt )
    {
        run(t,
            dt,
            &Simulation::compute_relaxed_increases,
            &Simulation::update_relaxed);
    }
    
    
    void save_profile( const string &fn, Real t ) const
    {
        ios::ocstream fp( fn, false);
        for(unit_t i=0; i <= imax; ++i )
        {
            fp("%g %g %g %g %g\n", X[i], h[i], w[i], -log10(h[i]), -log10( Kernel(t,X[i])));
        }
    }
    
    Real get_error() const throw()
    {
        Real err = 0;
        for( unit_t i=1; i < imax; ++i )
        {
            const Real Qw = h[i] * w[i];
            const Real tmp = Fabs( Qw - Kw) / Kw;
            if( tmp > err ) err = tmp;
        }
        return err;
    }
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Simulation);
};


static inline
double dt_round( double dt_max )
{
    const double dt_log = floor( log10(dt_max) );
    const double dt_one = floor(dt_max * pow(10.0,-dt_log));
    return dt_one * pow(10.0,dt_log);
}

#include "yocto/duration.hpp"
static
void process( eta &ETA, const double ratio, const char *kind , Real t)
{
    ETA(ratio);
    const string   percent = vformat("%7.1f%%", ratio*100 );
    const duration done(ETA.time_done);
    const duration left(ETA.time_left);
    const string   sdone = vformat("%02d:%02d:%02d",done.h,done.m,int(done.s));
    const string   sleft = vformat("%02d:%02d:%04.1f",left.h,left.m,left.s);
    
    std::cerr << kind << percent << " in " << sdone << " | ETA " << sleft << " @=t" << t << "  \r";
}

int main(int argc, char *argv[])
{
    const char *progname = _vfs::get_base_name( argv[0]);
    try
    {
        const Real alpha = 0.1;
        Simulation sim(500,1e-2);
        const Real dt_max   = alpha * (sim.dx*sim.dx) / max_of(sim.Dh,sim.Dw);
        const Real dt       = dt_round(dt_max);
        Real       dt_save  = 0.02;
        const Real t_run    = 5;
        Real       t        = 0;
        size_t     iter_max = 1+ceil(t_run/dt);
        size_t     every    = clamp<size_t>(1,dt_save/dt,iter_max);
        dt_save  = every * dt;
        while( 0 != (iter_max%every) ) ++iter_max;
        const size_t num_out = iter_max / every;
        
        std::cerr << "dt=" << dt << std::endl;
        std::cerr << "iter_max=" << iter_max << std::endl;
        std::cerr << "save every=" << every << ", dt_save=" << dt_save << std::endl;
        
        
        
        eta ETA;
        
        std::cerr << std::endl;
        std::cerr << "\t------------" << std::endl;
        std::cerr << "\tRelaxed" << std::endl;
        std::cerr << "\t------------" << std::endl;
        std::cerr << std::endl;
        
        sim.initialize();
        t=0;
        ETA.reset();
        sim.save_profile("h0.dat",0);
        for(size_t iter=1;iter<=iter_max;++iter)
        {
            sim.run_relaxed(t,dt);
            t = iter*dt;
            if( 0 == (iter%every) )
            {
                process(ETA,double(iter)/iter_max,"relaxed",t);
            }
        }
        std::cerr << std::endl;
        sim.save_profile("h1.dat",t);
        
        
        std::cerr << std::endl;
        std::cerr << "\t------------" << std::endl;
        std::cerr << "\tExplicit" << std::endl;
        std::cerr << "\t------------" << std::endl;
        std::cerr << std::endl;
        
        {
            sim.initialize();
            t=0;
            ETA.reset();
            for(size_t iter=1;iter<=iter_max;++iter)
            {
                sim.run_explicit(t,dt);
                t = iter*dt;
                if( 0 == (iter%every) )
                {
                    process(ETA,double(iter)/iter_max,"explicit",t);
                }
            }
            std::cerr << std::endl;
            sim.save_profile("h2.dat",t);
        }
        
        std::cerr << std::endl;
        std::cerr << "\t------------" << std::endl;
        std::cerr << "\tAll Done" << std::endl;
        std::cerr << "\t------------" << std::endl;
        std::cerr << std::endl;
        
        
        return 0;
    }
    catch( const exception &e )
    {
        std::cerr << "******** " << e.what() << std::endl;
        std::cerr << "******** " << e.when() << std::endl;
    }
    catch(...)
    {
        std::cerr << "******** Unhandled Exception in " << progname << "!" << std::endl;
    }
    
    return 1;
}
