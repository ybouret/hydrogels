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
        Initializer iniSide(L,"iniSide");
        std::cerr << iniSide << std::endl;
        
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
