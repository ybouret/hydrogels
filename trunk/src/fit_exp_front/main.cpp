#include "yocto/fs/vfs.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/exception.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/math/fit/lsf.hpp"
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
        const double slope = a[1];
        const double shift = a[2];
        return slope * (t-shift);
    }
    
    
    double compute2( double t, const array<double> &a)
    {
        const double D = a[1];
        return sqrt(D*max_of<double>(t,0));
    }
};

class Front
{
public:
    
    const size_t N;
    vector<double> t,x,x2,z2,z;
    vector<double> aorg;
    vector<double> aerr;
    const double   slope;
    const double   shift;
    
    explicit  Front( const string &fn ) :
    N(0),
    t(),
    x(),
    x2(),
    z2(),
    z(),
    aorg(2,0),
    aerr(2,0),
    slope(0),
    shift(0)
    {
        // load
        {
            ios::icstream fp(fn);
            data_set<double> ds;
            ds.use(3, t);
            ds.use(2, x);
            
            ds.load(fp,NULL,0,0);
            
            (size_t &)N = t.size();
            if(N<=2)
                throw exception("not enough points");
        }
        
        co_qsort(t, x);
        
        // finalize
        x2.make(N,0);
        z2.make(N,0);
        z.make(N,0);
        for(size_t i=1;i<=N;++i)
        {
            t[i] *= 10.0;    // in seconds
            x[i] *= 5.0/427; // in mm
            x2[i] = x[i] * x[i];
        }
        
        // prepare fit
        fit::sample<double> S( t, x2, z2);
        
        //-- eval slope
        aorg[1] = (x2[N]-x2[1])/(t[N]-t[1]);
        aorg[2] = 0;
        
        vector<bool> used(2,true);
        fit::lsf<double>    LeastSquares;
        LeastSquares.h = 1e-5;
        Diffusion                 diff;
        fit::lsf<double>::field   F( &diff, &Diffusion::compute );
        LeastSquares( S, F, aorg, used, aerr, NULL );
        if( S.status == fit:: success )
        {
            std::cerr << "Success: " << aorg << " +/-" << aerr << std::endl;
            
            
            (double&)slope  = aorg[1];
            (double&)shift  = aorg[2];
            const string output = fn + ".fit";
            const double slope_lo = slope - aerr[1];
            const double slope_hi = slope + aerr[1];
            
            for(size_t i=1;i<=N;++i)
            {
                z[i] = sqrt(max_of<double>(z2[i],0));
            }
            
            {
                const double Xerr = 5.0/427;
                ios::ocstream fp(output,false);
                fp("#t X Xerr F Ferr\n");
                fp("0 0 0 0 0\n");
                for(size_t i=1;i<=N;++i)
                {
                    const double tt = t[i] - shift;
                    const double zmax = sqrt( slope_hi * tt );
                    const double zmin = sqrt( slope_lo * tt );
                    const double zerr = (zmax-zmin)/2;
                    fp("%g %g %g %g %g\n", tt, x[i], Xerr, z[i], zerr);
                    
                }
            }
        }
        else
        {
            throw exception("Impossible fit!");
        }
    }
    
    
    typedef shared_ptr<Front> Ptr;
    
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Front);
};

int main( int argc, char *argv[] )
{
    const char *progname = vfs::get_base_name(argv[0]);
    try
    {
        if(argc<=1)
        {
            throw exception("usage: %s file1.txt ...",progname);
        }
        
        vector<Front::Ptr> fronts;
        
        for(int i=1;i<argc;++i)
        {
            Front::Ptr pFront( new Front( argv[i] ) );
            fronts.push_back(pFront);
        }
        
        
        vector<double> t;
        vector<double> x;
        vector<double> z;
        
        //t.push_back(0);
        //x.push_back(0);
        
        for(size_t i=1;i<=fronts.size();++i)
        {
            const Front &f = *fronts[i];
            for(size_t j=1;j<=f.N;++j)
            {
                t.push_back(f.t[j]-f.shift);
                x.push_back(f.x[j]);
            }
        }
        
        co_qsort(t, x);
        
        const size_t ntot = t.size();
        z.make(ntot,0);
        
        fit::sample<double> S(t,x,z);
        fit::lsf<double>    LeastSquares;
        LeastSquares.h = 1e-5;
        Diffusion                 diff;
        fit::lsf<double>::field   F( &diff, &Diffusion::compute2 );
        vector<double> aorg(1,0);
        vector<bool>   used(1,true);
        vector<double> aerr(1,0);
        
        aorg[1] = fronts[1]->slope;
        
        LeastSquares( S, F, aorg, used, aerr, NULL );
        
        if( S.status == fit:: success )
        {
            std::cerr << "aorg=" << aorg << " +/ " << aerr << std::endl;
            {
                ios::ocstream fp("fit.dat", false);
                fp("#t X F\n");
                fp("0 0 0\n");
                for(size_t i=1;i<=ntot;++i)
                {
                    fp("%g %g %g\n", t[i], x[i], z[i]);
                }
                
            }
        }
        
        
        
        
#if 0
        double slope     = 0;
        double slope_max = 0;
        double slope_min = 0;
        double error     = 0;
        size_t ntot  = 0;
        for(size_t i=1;i<=fronts.size();++i)
        {
            const Front &f = *fronts[i];
            slope += f.slope * f.N;
            ntot  += f.N;
            if(i==1)
            {
                slope_max = f.slope + f.aerr[1];
                slope_min = f.slope - f.aerr[1];
            }
            else
            {
                slope_max = max_of<double>(slope_max,f.slope+f.aerr[1]);
                slope_min = min_of<double>(slope_min,f.slope-f.aerr[1]);
            }
            error += f.N * (f.aerr[1]*f.aerr[1]);
        }
        slope /= ntot;
        error  = sqrt( error/ntot );
        
        std::cerr << "slope     = " << slope     << std::endl;
        std::cerr << "slope_max = " << slope_max << std::endl;
        std::cerr << "slope_min = " << slope_min << std::endl;
        std::cerr << "error     = " << error     << std::endl;
        
        vector<double> t(ntot,as_capacity);
        for(size_t i=1;i<=fronts.size();++i)
        {
            const Front &f = *fronts[i];
            const size_t n = f.N;
            for(size_t i=1;i<=n;++i) t.push_back(f.t[i]-f.aorg[2]);
        }
        quicksort(t);
        
        vector<double> t_err;
        vector<double> x_err;
        
        {
            ios::ocstream fp("fit.dat", false);
            fp("#t X n\n");
            fp("0 0 0\n");
            for(size_t i=1; i <= ntot; ++i )
            {
                const double zmax = sqrt(t[i]*slope_max);
                const double zmin = sqrt(t[i]*slope_min);
                const double zerr = (zmax-zmin)/2;
                fp("%g %g %g\n", t[i], sqrt( slope * t[i]), zerr);
                t_err.push_back(t[i]);
                x_err.push_back(zmin);
                
                t_err.push_back(t[i]);
                x_err.push_back(zmax);
            }
        }
#endif
        
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
