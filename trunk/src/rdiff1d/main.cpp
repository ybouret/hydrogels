#include "cell.hpp"

#include "yocto/fs/local-fs.hpp"
#include "yocto/exception.hpp"

#include "yocto/lua/lua-state.hpp"
#include "yocto/lua/lua-config.hpp"

static inline
bool is_curve( const vfs::entry &ep ) throw()
{
    return ep.has_extension("curve");
}

int  main(int argc, char *argv[] )
{
    const char *prog = vfs::get_base_name(argv[0]);
    try
    {
        if( argc <= 1 )
            throw exception("usage: %s file1.lua [file2.lua ...]",prog);
        
        vfs &fs = local_fs::instance();
        fs.create_sub_dir("data");
        fs.remove_files("data", is_curve);
        
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
        
        double dt = 0.1;
        double t  = 0.0;
        
        cell.init_all();
        cell.save_xy("data/v0.curve");
        
        for(int cycle=1;cycle<=200;++cycle)
        {
            cell.step(dt,t);
            t = cycle * dt;
            cell.save_xy("data/" + vformat("v%d.curve",cycle));
        }
        
        return 0;
    }
    catch(const exception &e)
    {
        std::cerr << "*** in " << prog  << std::endl;
        std::cerr << "*** "    << e.what() << std::endl;
        std::cerr << "*** "    << e.when() << std::endl;
    }
    catch(...)
    {
        std::cerr << "*** in " << prog << std::endl;
        std::cerr << "unhandled exception" << std::endl;
    }
    return 1;
}
