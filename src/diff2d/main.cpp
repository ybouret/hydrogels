#include "cell.hpp"

#include "yocto/lua/lua-state.hpp"
#include "yocto/lua/lua-config.hpp"

#include "yocto/string/vfs-utils.hpp"
#include "yocto/filesys/local-fs.hpp"

#include "yocto/ios/ocstream.hpp"
#include "yocto/code/utils.hpp"


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
        if( argc < 3 )
            throw exception("Usage: %s config.lua output_dir", progname);
        
        const string cfg    = argv[1];
        string       outdir = argv[2];
        Lua::State   VM;
        lua_State   *L = VM();
        Lua::Config::DoFile(L, cfg);
        
        //======================================================================
        //
        // Loading Simulation Parameters
        //
        //======================================================================
        Cell sim(L);

        const double dt = 0.1;
        sim.compute_fluxes();
        sim.compute_increases(0,dt);
        sim.vtk.save("increase.vtk", "increase", sim, sim.dvar, sim.vlayout);
        sim.update_fields(0.0);
        sim.vtk.save("v1.vtk","v1",sim,sim.var,sim.vlayout);
        return 0;
    }
    catch( const exception &e )
    {
        std::cerr << e.what() << std::endl;
        std::cerr << e.when() << std::endl;
    }
    catch(...)
    {
        std::cerr << "Unhandled Exception in " << progname << std::endl;
    }

    return 1;
    
}
