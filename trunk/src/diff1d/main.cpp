#include "initializer.hpp"
#include "cell.hpp"
#include "yocto/lua/lua-state.hpp"
#include "yocto/lua/lua-config.hpp"
#include "yocto/string/vfs-utils.hpp"

#include "yocto/ios/ocstream.hpp"

static inline
void save_h( const Cell &cell, const string &name )
{
    ios::ocstream fp( name, false );
    const Workspace &W = cell;
    const Array     &h = W["H+"].as<Array>();
    const Array     &X = cell.X;
    for( unit_t i=X.lower;i<=X.upper;++i)
    {
        fp("%g %g\n", X[i], h[i]);
    }
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
        Lua::State VM;
        lua_State *L = VM();
        Lua::Config::DoFile(L, cfg);

        
        //======================================================================
        //
        // Loading Simulation Parameters
        //
        //======================================================================
        Cell cell(L);
        cell.dt = 0.001;
        double shrink = -1;
        cell.initialize();
        save_h(cell, "h0.dat");
        
        cell.compute_fluxes();
        cell.compute_increases();
        cell.reduce();
        cell.find_shrink(shrink);
        cell.update();
        
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