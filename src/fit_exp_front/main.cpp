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

class Front
{
public:
    
    const string   filename;
    const size_t   N;
    vector<double> t,x,z;
    const double   Df;
    const double   t0;
    
    explicit  Front( const string &fn ) :
    filename(fn),
    N(0),
    t(),
    x(),
    z(),
    Df(0),
    t0(0)
    {
        // load
        {
            ios::icstream fp(filename);
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
        z.make(N,0);
        for(size_t i=1;i<=N;++i)
        {
            t[i] *= pix2tmx; // in seconds
            x[i] *= pix2pos; // in mm
        }
        
        // eval Df
        (double &)Df = (x[N]*x[N] - x[1]*x[1])/(t[N] - t[1]);
        
        
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
        
        LSF::samples samples;
        for(size_t i=1;i<=fronts.size();++i)
        {
            Front &f = *fronts[i];
            samples.append(f.t,f.x,f.z);
        }
        const size_t ns   = samples.size();
        const size_t npar = 2;
        const size_t nvar = 1 + ns;
        
        vector<double> a(nvar,0);
        a[nvar] = 0; // t0
        for(size_t i=1;i<=ns;++i)
        {
            LSF::sample &S =*samples[i];
            S.prepare(nvar, npar);
            matrix<double> &Gamma = S.Gamma;
            
            a[i] = fronts[i]->Df;
            Gamma[1][i]    = 1; //!< first param=Df for each
            Gamma[2][nvar] = 1; //!< last  param=t0 for all
        }
        Diffusion diff;
        LSF::function F( &diff, & Diffusion::compute);
        LSF lsf;
        lsf.h    = 1e-4;
        lsf.ftol = 1e-8;
        
        vector<double> aerr(nvar,0);
        vector<bool>   used(nvar,true);
        
        // first pass
        if( lsf(F,samples,a,used,aerr,NULL) != least_squares_failure )
        {
            std::cerr << "params=" << a << std::endl;
            std::cerr << "errors=" << aerr << std::endl;
            const double t_real = a[nvar];
            a[nvar]     = pix2tmx * floor( t_real/pix2tmx + 0.5);
            used[nvar] = false;
        }
        else
        {
            throw exception("Couldn't fit!!!");
        }
        std::cerr << "Exp t0=" << a[nvar] << std::endl;
        
        
        
        const string ext = vformat(".fit%u", unsigned(ns) );
        
        if( lsf(F,samples,a,used,aerr,NULL) != least_squares_failure )
        {
            std::cerr << "params=" << a << std::endl;
            std::cerr << "errors=" << aerr << std::endl;
            
            double Dmin = 0;
            double Dmax = 0;
            double Emax = 0;
            vector<double> t;
            for(size_t i=1;i<=ns;++i)
            {
                Front &f = *fronts[i];
                (double &)(f.Df) = a[i];
                (double &)(f.t0) = a[nvar];
                if(i==1)
                {
                    Dmin = f.Df;
                    Dmax = f.Df;
                }
                else
                {
                    Dmin = min_of(Dmin,f.Df);
                    Dmax = max_of(Dmax,f.Df);
                }
                Emax = max_of(Emax,Fabs(aerr[i]));
                
                t.reserve(f.N);
                for(size_t j=1;j<=f.N;++j)
                {
                    t.push_back(f.t[j]);
                }
                
                {
                    ios::ocstream fp( f.filename + ext, false);
                    fp("#t pos fit err\n");
                    //fp("0 0 0 0\n");
                    for(size_t j=1;j<=f.N;++j)
                    {
                        const double tj = f.t[j];
                        const double xj = f.x[j];
                        const double dt = tj - f.t0;
                        const double Fj = sqrt( f.Df * dt );
                        const double Fmax  = sqrt( (f.Df+aerr[i]) * dt );
                        const double Fmin = sqrt( (f.Df-aerr[i]) * dt );
                        const double Ferr = Fmax-Fmin;
                        fp("%g %g %g %g\n",dt,xj,Fj,Ferr);
                    }
                }
                {
                    ios::ocstream fp( f.filename + ext + ".pix", false);
                    for(size_t j=1;j<=f.N;++j)
                    {
                        const double Yj = (f.t[j] / pix2tmx);
                        const double dt = f.t[j] - f.t0;
                        const double Xj = sqrt( f.Df * dt )/pix2pos;
                        fp("%g %g %g\n", f.x[j]/pix2pos, Yj, Xj);
                    }
                }
                
            }
            std::cerr << "Dmin=" << Dmin << std::endl;
            std::cerr << "Dmax=" << Dmax << std::endl;
            const double Df   = 0.5 * (Dmax+Dmin);
            const double Derr = 0.5 * (Dmax-Dmin) + Emax;
            std::cerr << "Df  =" << Df << std::endl;
            std::cerr << "Derr=" << Derr  << std::endl;
            
            quicksort(t);
            
            {
                ios::ocstream fp("water-front" + ext,false);
                fp("#t pos err\n");
                for(size_t j=1;j<=t.size();++j)
                {
                    const double dt = t[j]-a[nvar];
                    const double xx = sqrt( Df * dt );
                    const double xxmax = sqrt( (Df+Derr) * dt );
                    const double xxmin = sqrt( (Df-Derr) * dt );
                    const double xerr  = (xxmax-xxmin)/2;
                    fp("%g %g %g\n", dt, xx, xerr);
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
