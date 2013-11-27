#include "cell.hpp"

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
            throw exception("usage: %s file1.lua [file2.lua ...]",prog);
        
        //______________________________________________________________________
        //
        // Loading Lua resources
        //______________________________________________________________________

        Lua::State VM;
        lua_State *L = VM();
        for(int i=1;i<argc;++i)
        {
            std::cerr << "... loading file '" << argv[1] << "'" << std::endl;
            Lua::Config::DoFile(L, argv[i]);
        }
        
        //______________________________________________________________________
        //
        // Massive Initialization
        //______________________________________________________________________
        
        Cell cell(L);
        
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
