#include "cell.hpp"

#include "yocto/lua/lua-state.hpp"
#include "yocto/lua/lua-config.hpp"

#include "yocto/string/vfs-utils.hpp"
#include "yocto/filesys/local-fs.hpp"

#include "yocto/ios/ocstream.hpp"
#include "yocto/code/utils.hpp"
#include "yocto/filesys/local-fs.hpp"
#include "yocto/math/round.hpp"
#include "yocto/code/rand.hpp"
#include "yocto/eta.hpp"

using namespace filesys;


static inline bool is_data_file( const vfs::entry &ep )
{
    return ep.is_reg && ep.extension != 0 && 0 == strcmp( ep.extension, "vtk" );
}

void save_data( const string &outdir, const Cell &sim, unsigned cycle )
{
    const string filename = outdir + vformat( "f%u.vtk", cycle );
    sim.vtk.save(filename,"v1",sim,sim.var,sim.vlayout);
    
}

int main( int argc, char *argv[] )
{
    const char *progname = _vfs::get_base_name( argv[0]);
    try
    {
        _rand.seed( wtime::seed() );
        
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
        
        for( unit_t j=sim.lower.y; j <= sim.upper.y ; ++j)
        {
            for( unit_t i=sim.lower.x; i <= sim.upper.x; ++i )
            {
                sim.porosity[j][i] = 0.5 + 0.5 * alea<double>();
            }
        }
        
        //======================================================================
        //
        // initialize output database
        //
        //======================================================================
        vfs & fs = local_fs::instance();
        _vfs::as_directory(outdir);
        fs.create_sub_dir(outdir);
        std::cerr << "cleaning " << outdir << std::endl;
        { vfs::entry::callback filter( cfunctor( is_data_file ) ); fs.remove_files(outdir, filter); }
        
        //======================================================================
        //
        // Time parameters
        //
        //======================================================================
        const double t_run   = Lua::Config::Get<lua_Number>(L,"t_run");
        const double alpha   = Lua::Config::Get<lua_Number>(L,"alpha");
        const double dt_save = Lua::Config::Get<lua_Number>(L,"dt_save");
        double Dmax = 0;
        {
            const library &l = sim;
            for( library::const_iterator i=l.begin(); i != l.end(); ++i )
            {
                const species     &sp = **i;
                const SpeciesData &sd = sp.get<SpeciesData>();
                const double D = sd.D;
                if( D > Dmax ) Dmax = D;
            }
        }
        if( Dmax <= 0 ) throw exception("Invalid Diffusion Coefficients");
        std::cerr << "Dmax   = " << Dmax << std::endl;
        const double dt_raw = alpha * min_of( sim.del_sq.x, sim.del_sq.y ) / Dmax;
        std::cerr << "dt_raw = " << dt_raw << std::endl;
        const double dt = log_round(dt_raw);
        std::cerr << "dt     = " << dt << std::endl;
        size_t         iter    = 1 + ceil( t_run/dt );
        const size_t   every   = clamp<size_t>(1,floor(dt_save/dt),iter);
        while( 0 != (iter%every) ) ++iter;
        std::cerr << "#steps = " << iter << std::endl;
        
        size_t isave = 0;
        save_data(outdir, sim,isave);
        for(unsigned cycle=1;cycle<=iter;++cycle)
        {
            const double t = (cycle-1) * dt;
            sim.compute_increases(t,dt);
            sim.update_fields(t);
            if( 0 == ( cycle % every ) )
            {
                ++isave;
                save_data(outdir, sim, isave);
                std::cerr << "t=" << t << " / " << t_run << std::endl;
            }
        }
        
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
