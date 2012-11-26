#include "yocto/string/vfs-utils.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/exception.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/math/fit/lsf.hpp"
#include "yocto/code/utils.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/stat/descr.hpp"

#include <iostream>

using namespace yocto;
using namespace math;


struct diffusion
{
    double compute( double x, const array<double> &a )
    {
        return sqrt( a[1] * max_of<double>(0,x) );
    }
    
    bool cb(const fit::sample<double>     &s, 
            const fit::lsf<double>::field &f,
            const array<double>           &a )
	{
        std::cerr << std::endl;
        std::cerr << "\tDistance:   " << s.D << std::endl;
        std::cerr << "\tVariables = " << a << std::endl;
        return true;
    }
    
};


int main( int argc, char *argv[] )
{
    const char *progname = _vfs::get_base_name(argv[0]);
    try
    {
        
        if(argc<=3)
            throw exception("usage: %s data_file factor cut", progname);
        const string filename = argv[1];
        const double factor   = strconv::to_double( argv[2], "factor" );
        const double cut      = strconv::to_double( argv[3], "cut" );
        vector<double> t;
        vector<double> y;
        std::cerr << "-- Loading Data" << std::endl;
        {
            data_set<double> ds;
            ds.use(1,t);
            ds.use(2,y);
            ios::icstream fp( filename );
            ds.load(fp);
            for( size_t i=y.size();i>0;--i) 
                y[i] *= factor;
            if( cut > 0 )
            {
                while( (y.size() > 0) && (y.back() > cut) )
                {
                    y.pop_back();
                    t.pop_back();
                }
            }
        }
        
        const size_t n = y.size();
        std::cerr << "-- Processing #data= " << n << std::endl;
        
        
        vector<double> z(n,0.0);
        
        std::cerr << "-- Prepare the Fit Function" << std::endl;
        diffusion                   diff;
        fit::lsf<double>::field     F( &diff, &diffusion::compute );
        fit::lsf<double>::callback  G( &diff, &diffusion::cb      );
        
        std::cerr << "-- Prepare the Fit Variables" << std::endl;
        const size_t nv = 1;
        const double D0 = 1e-9;
        vector<double> aorg(nv,1e-9);
        vector<bool>   used(nv,true);
        vector<double> aerr(nv,-1);
        
        std::cerr << "-- Prepare the Fit Sample" << std::endl;
        fit::sample<double> S( t, y, z );
        
        std::cerr << "-- Prepare the LeastSquares" << std::endl;
        fit::lsf<double>    LeastSquares;
        LeastSquares.h = D0 * 1e-3;
        LeastSquares( S, F, aorg, used, aerr, &G );
        if( S.status == fit:: success )
        {
            std::cerr << std::endl;
            std::cerr << "\tD    = " << aorg[1] << " +/- " << aerr[1]/2 << " m^2/s" << std::endl;
            const double R = compute_correlation(y, z);
            std::cerr << "\tcorr = " << R << std::endl;
            std::cerr << "-- Saving Data" << std::endl;
            {
                const string outfile = filename + ".fit";
                ios::ocstream fp( outfile, false );
                for( size_t i=1; i <= n; ++i )
                {
                    fp("%.15g %.15g %.15g\n", t[i], y[i], z[i] );
                }
            }
            std::cerr << "-- Saving Log" << std::endl;
            {
                const string outfile = filename + ".log";
                ios::ocstream fp( outfile, false );
                fp("D      = %.15g +/- %.15g m^2/s\n", aorg[1], aerr[1]/2);
                fp("corr   = %.8f\n", R );
                fp("factor = %.15g\n", factor);
                fp("Cut at = %.15g m\n", cut);
            }
            
        }
        else 
        {
            std::cerr << "Failure" << std::endl;
        }
        return 0;
    }
    catch( const exception &e )
    {
        std::cerr << "in " << progname << std::endl;
        std::cerr << e.what() << std::endl;
        std::cerr << e.when() << std::endl;
    }
    catch(...)
    {
        std::cerr << "unhandled exception in " << progname << std::endl;
    }
	return 1;
    
    
}