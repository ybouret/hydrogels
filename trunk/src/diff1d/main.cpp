#include "initializer.hpp"
#include "cell.hpp"
#include "yocto/lua/lua-state.hpp"
#include "yocto/lua/lua-config.hpp"
#include "yocto/string/vfs-utils.hpp"

#include "yocto/ios/ocstream.hpp"
#include "yocto/code/utils.hpp"

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
void save_front( const double pH, const Cell &cell, const double t)
{
    ios::ocstream fp( "front.dat", t>0 );
    fp("%g %g\n", t, cell.locate(pH) );
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
        const double pH_front = 6;
        
        //======================================================================
        //
        // Adjusting iterations
        //
        //======================================================================
        const double t_run = Lua::Config::Get<double>(L,"t_run");
        const size_t iter  = 1 + size_t(ceil(t_run/dt));
        const double dt_save = Lua::Config::Get<double>(L,"dt_save");
        const size_t every = clamp<size_t>(1,floor(dt_save/dt),iter);
        std::cerr << "###  dt_max=" << dt_max << std::endl;
        std::cerr << "###      dt=" << dt << std::endl;
        std::cerr << "###   t_run=" << t_run << std::endl;
        std::cerr << "###    iter=" << iter << std::endl;
        std::cerr << "### dt_save=" << dt_save << std::endl;
        std::cerr << "###   every=" << every << " / dt_out=" << every * dt << std::endl;
        double t=0;
        cell.initialize();
        
        save_h(cell, "h0.dat");
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
                save_front(pH_front, cell, t);
            }
            std::cerr << "t= " << t << "\r"; std::cerr.flush();
        }
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