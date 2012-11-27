#include "initializer.hpp"
#include "yocto/lua/lua-state.hpp"
#include "yocto/lua/lua-config.hpp"
#include "yocto/string/vfs-utils.hpp"
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
        Library     lib(L);
        ChemSys     cs(lib,L);
        Initializer ini("ini_left",lib,L);
        ini(cs,0.0);
        
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