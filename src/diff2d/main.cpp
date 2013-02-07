#include "chemsys.hpp"
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
        Library lib(L);
        ChemSys cs(lib,L);
        
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
