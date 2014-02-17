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
    double Compute( double t, const array<double> &a )
    {
        const double Df = a[2];
        const double t0 = a[1];
        return Df * (t-t0);
    }
};

static inline double __sqr(double x) throw() { return x*x; }


int main( int argc, char *argv[] )
{
    const char *progname = vfs::get_base_name(argv[0]);
    try
    {
        if(argc<=2)
        {
            throw exception("usage: %s file.txt t_mini", progname);
        }
        
        const string filename = argv[1];
        const double t_mini   = strconv::to<double>(argv[2],"t_mini");
        
        vector<double> t;
        vector<double> x;
        {
            data_set<double> ds;
            vector<double>   tt;
            vector<double>   xx;
            ds.use(1,tt);
            ds.use(2,xx);
            ios::icstream fp(filename);
            ds.load(fp);
            const size_t nt = tt.size();
            t.reserve(nt);
            x.reserve(nt);
            
            for(size_t i=1;i<=nt;++i)
            {
                if( xx[i] >= 0 && tt[i] >= t_mini)
                {
                    t.push_back(tt[i]);
                    x.push_back(xx[i]);
                }
            }
        }
        
        co_qsort(t, x);
        const size_t n    = t.size();
        const size_t nvar = 2;
        
        vector<double> aorg(nvar,0.0);
        vector<double> aerr(nvar,0.0);
        vector<bool>   used(nvar,true);
        vector<double> x2(n,0.0);
        vector<double> z(n,0.0);
        
        least_squares<double>::samples samples;
        samples.append(t,x2,z);
        for(size_t i=1;i<=n;++i)
            x2[i] = __sqr(x[i]);
        
        least_squares<double>::sample &sample = *samples[1];
        
        double &Df = aorg[2];
        double &t0 = aorg[1];
        
        double &Derr = aerr[2];
        double &terr = aerr[1];
        
        sample.polynomial(aorg,used,aerr);
        {
            const double inter = aorg[1];
            const double slope = aorg[2];
            Df = slope;
            t0 = -inter/slope;
            std::cerr << "Df=" << Df << ", t0=" << t0 << std::endl;;
            
            ios::ocstream fp("output.dat", false);
            for(size_t i=1;i<=n;++i)
            {
                const double dt = t[i] - t0;
                fp("%g %g %g\n", t[i], x[i], sqrt(Df*dt) );
            }
        }
        
        least_squares<double> LSF;
        LSF.h = 1e-5;
        Diffusion diff;
        least_squares<double>::function F(&diff, &Diffusion::Compute);
        if( LSF(F,samples,aorg,used,aerr) != least_squares_failure )
        {
            std::cerr << "Df=" << Df << " +/-" << Derr << std::endl;
            std::cerr << "t0=" << t0 << " +/-" << terr << std::endl;
        
            {
                const string  output = filename + ".fit";
                ios::ocstream fp(output,false);
                fp("#t x fit\n");
                for(size_t i=1;i<=n;++i)
                {
                    const double dt = t[i] - t0;
                    z[i] = sqrt(max_of<double>(Df*dt,0));
                    if(dt>=0)
                    {
                        fp("%g %g %g\n", dt, x[i], z[i]);
                    }
                }
                
            }
            
            {
                const string output = filename + ".log";
                ios::ocstream fp(output,false);
                fp("Df= %.8g +/- %.8g\n", Df, Derr);
                fp("t0= %.8g +/- %.8g\n", t0, terr);
                fp("corr= %.8g\n", compute_correlation(x, z));
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
