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
        
        const double alpha = 0.1;
        const double dx_min = cell.dx_min();
        const double D_max = cell.D_max();
        if(D_max<=0)
            throw exception("No Diffusion!");
        
        const double dt_max = alpha * dx_min * dx_min/D_max;
        std::cerr << "dt_max=" << dt_max << std::endl;
        
        
        
        double dt = dt_round(dt_max);
        std::cerr << "dt=" << dt << std::endl;
        double t  = 0.0;
        
        cell.init_all();
        cell.save_xy("data/v0.curve");
        
        for(int cycle=1;t<=10;++cycle)
        {
            cell.step(dt,t);
            t = cycle * dt;
            cell.save_xy("data/" + vformat("v%d.curve",cycle));
            std::cerr << "t=" << t << std::endl;
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
