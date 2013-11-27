#include "initializer.hpp"
#include "yocto/fs/vfs.hpp"
#include "yocto/exception.hpp"

#include "yocto/lua/lua-state.hpp"
#include "yocto/lua/lua-config.hpp"

int  main(int argc, char *argv[] )
{
    const char *prog = vfs::get_base_name(argv[0]);
    try
    {
        if( argc <= 1 )
            throw exception("usage: %s config.lua",prog);
        
        std::cerr << "... loading file '" << argv[1] << "'" << std::endl;
        Lua::State VM;
        lua_State *L = VM();
        Lua::Config::DoFile(L, argv[1]);
        
        Collection  lib(L);
        
        ChemSys     cs(L,lib);
        std::cerr << cs << std::endl;
        
        Initializer ini_side(L,"ini_side",lib);
        std::cerr << "ini@side=" << std::endl << ini_side << std::endl;
       
        Initializer ini_core(L,"ini_core",lib);
        std::cerr << "ini@core=" << std::endl << ini_core << std::endl;
        
        
        ini_side(cs,lib,0.0);
        
        chemical::solution S(lib);
        S.load(cs.C);
        std::cerr << "side=" << S << std::endl;
        std::cerr << "pH=" << S.pH() << std::endl;
        return 0;
    }
    catch(const exception &e)
    {
        std::cerr << "*** in " << prog  << std::endl;
        std::cerr << "*** " << e.what() << std::endl;
        std::cerr << "*** " << e.when() << std::endl;
    }
    catch(...)
    {
        std::cerr << "*** in " << prog << std::endl;
        std::cerr << "unhandled exception" << std::endl;
    }
    return 1;
}
