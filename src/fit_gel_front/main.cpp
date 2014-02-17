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

static inline double __sqr(double x) throw() { return x*x; }


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
            
            co_qsort(t,x);
            const size_t   n = t.size();
            vector<double> z(n,0);
            
            least_squares<double>::samples samples;
            samples.append(t,x,z);
            
            least_squares<double> LSF;
            
            // variables
            vector<double> a(2,0.0);
            vector<bool>   used(2,true);
            vector<double> aerr(2,0.0);
            
            //! diffusion coefficient estimation
            double &Df = a[1];
            double &t0 = a[2];
            
            {
                least_squares<double>::sample &sample = *samples[1];
                numeric<double>::function      transf = cfunctor(__sqr);
                sample.polynomial(a, used, aerr, &transf);
                const double slope = a[2];
                const double inter = a[1];
                std::cerr << "slope=" << slope << ", inter=" << inter << std::endl;
                Df = slope;
                t0 = -inter/Df;
                std::cerr << "Df=" << Df << ", t0=" << t0 << std::endl;
            }
            t0 = 0;
        
            Diffusion diff;
            least_squares<double>::function F( &diff, & Diffusion::compute );
          
            if(LSF(F,samples,a,used,aerr,NULL) != least_squares_failure )
            {
                std::cerr << "Df=" << Df << " +/-" << aerr[1] << std::endl;
                std::cerr << "t0=" << t0 << " +/-" << aerr[2] << std::endl;
                ios::ocstream fp("output.dat",false);
                for(size_t i=1;i<=n;++i)
                {
                    fp("%g %g %g\n", t[i], x[i], (z[i]));
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
