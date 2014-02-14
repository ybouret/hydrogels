#include "yocto/fs/vfs.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/exception.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/math/fit/least-squares.hpp"
#include "yocto/code/utils.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/stat/descr.hpp"
#include "yocto/sort/quick.hpp"
#include "yocto/math/fcn/functions.hpp"
#include "yocto/ptr/shared.hpp"

#include <iostream>

using namespace yocto;
using namespace math;


struct  Diffusion
{
    double compute( double t, const array<double> &a )
    {
        const double Df = a[1];
        const double t0 = a[2];
        return sqrt( max_of<double>(Df * (t-t0),0 ) );
    }
};

typedef least_squares<double> LSF;

static const double pix2tmx = 10.0;
static const double pix2pos = 5.0/427;


int main( int argc, char *argv[] )
{
    const char *progname = vfs::get_base_name(argv[0]);
    try
    {
        if(argc<=1)
        {
            throw exception("usage: %s file1.txt ...",progname);
        }
        
        for(int na=1;na<argc;++na)
        {
            const string filename = argv[na];
            vector<double> t;
            vector<double> x;
            {
                data_set<double> ds;
                vector<double> tt;
                vector<double> xx;
                ds.use(1,tt);
                ds.use(2,xx);
                ios::icstream fp(filename);
                ds.load(fp);
                for(size_t i=1;i<=tt.size();++i)
                {
                    if( xx[i] >= 0 )
                    {
                        t.push_back(tt[i]);
                        x.push_back(xx[i]);
                        
                    }
                }
            }
            const size_t n = t.size();
            vector<double> x2(n,0);
            
            {
                ios::ocstream fp("output.dat",false);
                for(size_t i=1;i<=n;++i)
                {
                    fp("%g %g\n", t[i], x[i]);
                }
            }
        }
        
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
