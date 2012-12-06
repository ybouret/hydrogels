#include "initializer.hpp"
#include "cell.hpp"
#include "yocto/lua/lua-state.hpp"
#include "yocto/lua/lua-config.hpp"
#include "yocto/string/vfs-utils.hpp"

#include "yocto/ios/ocstream.hpp"
#include "yocto/code/utils.hpp"

#include "yocto/math/fit/lsf.hpp"

static inline
void save_h( const Cell &cell, const string &name )
{
    ios::ocstream fp( name, false );
    const Workspace &W = cell;
    const Array     &h = W["H+"].as<Array>();
    const Array     &X = cell.X;
    for( unit_t i=X.lower;i<=X.upper;++i)
    {
        //fp("%g %g\n", X[i], (h[i]));
        fp("%g %g\n", X[i], -log10(h[i]));
    }
}

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

static inline
double save_front( const double pH, const Cell &cell, const double t)
{
    const string fn = vformat("front%d.dat", int( pH*100) );
    ios::ocstream fp( fn, t>0 );
    const double x = cell.locate(pH);
    fp("%g %g\n", t, x );
    return x;
}

double FitFunction( double t, const array<double> &a )
{
    return a[1] * t;
}

int main( int argc, char *argv[] )
{
    const char *progname = _vfs::get_base_name( argv[0]);
    try
    {
        //======================================================================
        //
        // Loading Lua Parameters
        //
        //======================================================================
        if( argc < 2 )
            throw exception("Usage: %s config.lua", progname);
        
        const string cfg = argv[1];
        Lua::State   VM;
        lua_State   *L = VM();
        Lua::Config::DoFile(L, cfg);
        
        
        //======================================================================
        //
        // Loading Simulation Parameters
        //
        //======================================================================
        Cell cell(L);
        
        //======================================================================
        //
        // Adjusting time step
        //
        //======================================================================
        const double dt_max = cell.max_dt();
        double dt=dt_round(dt_max);
        
        //======================================================================
        //
        // Look up front ?
        //
        //======================================================================
        double pH_front    = 7;
        bool   build_front = false;
        {
            lua_getglobal(L, "pH_front");
            if( lua_isnumber(L, -1))
            {
                build_front = true;
                pH_front    = lua_tonumber(L, -1);
                std::cerr << "Building Front @pH=" << pH_front << std::endl;
            }
        }
        
        //======================================================================
        //
        // Fitting
        //
        //======================================================================
        fit::lsf<double>    Fit;
        vector<double>      fX,fY,fZ;
        fit::sample<double> Sample(fX,fY,fZ);
        vector<double>      coef;
        vector<bool>        used;
        vector<double>      aerr;
        
        const bool          fit = Lua::Config::Get<bool>(L,"fit");
        std::cerr << ( fit ? "Fitting Enabled" : "Fitting Disabled") << std::endl;
        fit::lsf<double>::field func = cfunctor2(FitFunction);
        
        //======================================================================
        //
        // Adjusting iterations
        //
        //======================================================================
        const double t_run = Lua::Config::Get<double>(L,"t_run");
        size_t       iter  = 1 + size_t(ceil(t_run/dt));
        const double dt_save = Lua::Config::Get<double>(L,"dt_save");
        const size_t every = clamp<size_t>(1,floor(dt_save/dt),iter);
        while( 0 != (iter%every) ) ++iter;
        const size_t nsave = iter/every;
        std::cerr << "###  dt_max=" << dt_max << std::endl;
        std::cerr << "###      dt=" << dt << std::endl;
        std::cerr << "###   t_run=" << t_run << std::endl;
        std::cerr << "###    iter=" << iter << std::endl;
        std::cerr << "###   t_run=" << iter * dt << std::endl;
        std::cerr << "### dt_save=" << dt_save << std::endl;
        std::cerr << "###   every=" << every << " / dt_out=" << every * dt << std::endl;
        std::cerr << "###   nsave=" << nsave << std::endl;
        double t=0;
        cell.initialize();
        
        if( fit )
        {
            fX.reserve(nsave);
            fY.reserve(nsave);
            fZ.reserve(nsave);
            coef.make(1,0);
            used.make(1,true);
            aerr.make(1,0);
        }
        
        
        save_h(cell, "h0.dat");
        if( build_front )
            save_front( pH_front, cell, t );
        
        
        double ell = 0;
        size_t nst = 0;
        for( size_t i=1; i <= iter; ++i )
        {
            ell += cell.step(t,dt);
            ++nst;
            t = i * dt;
            if( 0 == (i%every) )
            {
                std::cerr << "t= " << t << "\r"; std::cerr.flush();
                if(build_front)
                {
                    const double x = save_front(pH_front, cell, t);
                    if(fit)
                    {
                        fX.push_back(t);
                        fY.push_back(x*x);
                        fZ.push_back(0);
                        Fit(Sample,func,coef,used,aerr);
                        std::cerr << std::endl << "D=" << coef[1] * 1e8<< ", err=" << aerr[1]*1e8 << std::endl;
                    }
                }
            }
        }
        std::cerr << std::endl;
        std::cerr << "<steps/s>=" << nst/ell << std::endl;
        
        save_h(cell,"h1.dat");
        return 0;
    }
    catch( const exception &e )
    {
        std::cerr << e.what() << std::endl;
        std::cerr << e.when() << std::endl;
    }
    catch(...)
    {
        std::cerr << "Unhandled Exception!" << std::endl;
    }
    
    return 1;
}