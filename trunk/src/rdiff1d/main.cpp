#include "cell.hpp"

#include "yocto/fs/local-fs.hpp"
#include "yocto/exception.hpp"

#include "yocto/lua/lua-state.hpp"
#include "yocto/lua/lua-config.hpp"

#include "yocto/code/utils.hpp"
#include "yocto/eta.hpp"
#include "yocto/duration.hpp"

#include "yocto/ios/ocstream.hpp"


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
    const char *progname = vfs::get_base_name(argv[0]);
    try
    {
        if( argc <= 1 )
            throw exception("usage: %s file1.lua [file2.lua| -f \"code\"]",progname);
        
        vfs &fs = local_fs::instance();
        fs.create_sub_dir("data");
        fs.remove_files_with_extensions("data", "curve");
        
        //______________________________________________________________________
        //
        // Loading Lua resources
        //______________________________________________________________________
        
        Lua::State VM;
        lua_State *L = VM();
        for(int i=1;i<argc;++i)
        {
            const string args = argv[i];
            if( args == "-f" )
            {
                if(++i>=argc)
                    throw exception("missing code after flag!!!");
                const string code = argv[i];
                std::cerr << "... compiling code '" << code << "'" << std::endl;
                Lua::Config::DoString(L, code);
            }
            else
            {
                std::cerr << "... loading file '" << args << "'" << std::endl;
                Lua::Config::DoFile(L, argv[i]);
            }
            
        }
        
        //______________________________________________________________________
        //
        // Massive Initialization
        //______________________________________________________________________
        
        Cell cell(L);
        
        
        if(cell.search_front)
        {
            ios::ocstream::overwrite(cell.search_output);
        }
        else
        {
            try { fs.remove_file(cell.search_output); } catch(...) {}
        }
        //______________________________________________________________________
        //
        // time control parameters
        //______________________________________________________________________
        const double alpha  = clamp<double>(0,Lua::Config::Get<lua_Number>(L, "alpha"),0.5);
        const double dx_min = cell.dx_min();
        const double D_max  = cell.D_max();
        if(D_max<=0)
            throw exception("No Diffusion!");
        
        const double Tmax = Lua::Config::Get<lua_Number>(L, "Tmax");
        
        
        const double dt_alpha = dt_round(alpha * dx_min * dx_min/D_max);
        std::cerr << "dt_alpha=" << dt_alpha << std::endl;
        
        const double dt_user = Lua::Config::Get<lua_Number>(L, "dt" );
        std::cerr << "dt_user=" << dt_user << std::endl;
        const double dt = min_of(dt_user,dt_alpha);
        
        std::cerr << "dt=" << dt << std::endl;
        
        const double   dt_sav = max_of<double>(0,Lua::Config::Get<lua_Number>(L, "save"));
        const unsigned every  = max_of<unsigned>(1,ceil(dt_sav/dt));
        std::cerr << "saving every " << every << " steps= " << dt * every << " / " << dt_sav << std::endl;
        //______________________________________________________________________
        //
        // Let us go
        //______________________________________________________________________
        eta    prog;
        double t  = 0.0;
        
        std::cerr.flush();
        prog.reset();
        cell.init_all();
        int cycle = 0;
        int isave = 0;
        double old_tmx = prog.now();
        for(;t<=Tmax;++cycle)
        {
            t=cycle*dt;
            if( 0 ==( cycle % every) )
            {
                cell.save_xy("data/" + vformat("v%d.curve",isave++));
                prog(t/Tmax);
                const duration tmx_done( prog.time_done );
                const duration tmx_left( prog.time_left );
                const double   new_tmx = prog.now();
                const double   fps = every / (new_tmx - old_tmx);
                old_tmx = new_tmx;
                fprintf(stderr,"Time: %02uD%02uH%02uM%04.1fs | ETA:  %02uD%02uH%02uM%04.1fs | sim= %.1fs @ %.1f fps  \r",
                        tmx_done.d, tmx_done.h, tmx_done.m, tmx_done.s,
                        tmx_left.d, tmx_left.h, tmx_left.m, tmx_left.s,
                        t,
                        fps
                        );
                
                fflush(stderr);
                if(cell.search_front)
                {
                    if(cell.search_value<=0)
                    {
                        double xx=0;
                        double yy=0;
                        if( cell.find_inflection(xx,yy) )
                        {
                            ios::ocstream fp(cell.search_output,true);
                            fp("%g %.15e %.15e\n",t,xx,yy);
                        }
                        
                    }
                    else
                    {
                        double pos = 0;
                        if( cell.find_front(pos) )
                        {
                            ios::ocstream fp(cell.search_output,true);
                            fp("%g %.15e\n",t,pos);
                        }
                    }
                }
            }
            cell.step(dt,t);
        }
        std::cerr << std::endl;
        return 0;
    }
    catch(const exception &e)
    {
        std::cerr << "*** in " << progname  << std::endl;
        std::cerr << "*** "    << e.what() << std::endl;
        std::cerr << "*** "    << e.when() << std::endl;
    }
    catch(...)
    {
        std::cerr << "*** in " << progname << std::endl;
        std::cerr << "unhandled exception" << std::endl;
    }
    return 1;
}
