#include "yocto/math/fcn/zfind.hpp"
#include "yocto/math/fcn/functions.hpp"

#include "yocto/fs/vfs.hpp"
#include "yocto/exception.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/sequence/vector.hpp"

using namespace yocto;

static const double Kw  = 1e-14;
static double       pH0 = 2;
static double       h0  = pow(10,-pH0);
static double       w0  = Kw/h0;
static double       s0  = h0+w0;
static double       d0  = h0-w0;

static double       pH1 = 10;
static double       h1  = pow(10,-pH1);
static double       w1  = Kw/h1;
static double       s1  = h1+w1;
static double       d1  = h1-w1;

static double D = 5e-9;

static inline double get_arg(double t,double x )
{
    return x/sqrt(4*D*t);
}

static inline double d(double t, double x)
{
    if(x<=0)
        return d0;
    else
    {
		const double arg = get_arg(t,x);
        return d1 + (d0-d1) * math::qerfc(arg);
	}
}

static inline double s(double t, double x)
{
    if(x<=0)
        return s0;
    else
    {
		const double dtx =  d(t,x);
        return sqrt(s1*s1 + dtx*dtx - d1*d1);
    }
}

static inline double get_h(double t, double x )
{
    return 0.5 * ( d(t,x) + s(t,x) );
}

static inline
double get_alpha( double h )
{
    const double dd  = h - Kw/h;
    const double rho = (dd - d1)/(d0-d1);
    if( rho <=0 || rho >= 1)
        throw exception("Invalid h=%g", h);
    //std::cerr << "h=" << h << ", rho=" << rho << std::endl;
    const double q   = math::iqerf(1.0-rho);
    return 4*q*q;
}

int main(int argc, char *argv[])
{
    const char *progname = vfs::get_base_name(argv[0]);
    try
    {
        
        vector<double> t;
        t.push_back(25);
        t.push_back(50);
        t.push_back(100);
        t.push_back(200);
        t.push_back(400);
        t.push_back(800);
        t.push_back(1600);

        
        double dx = 0.0001;
        {
			ios::ocstream fp( "profils.dat", false );
            fp("0 ");
            for(size_t i=1; i<=t.size();++i)
            {
                fp(" %g", pH0);
            }
            fp("\n");
            
	        for( double x=dx; x <= 0.02; x += dx )
	        {
                fp("%g",x);
                for(size_t i=1; i<=t.size();++i)
                {
                    fp(" %g", -log10(get_h(t[i],x)) );
                }
                fp("\n");
			}
		}
        
        {
            ios::ocstream fp( "Deff.dat", false);
            const double dph = 0.01;
            for( double pH = pH0+dph; pH < pH1-dph; pH += dph)
            {
                const double alpha = get_alpha( pow(10.0,-pH) );
                fp("%g %g\n", pH, alpha );
            }
        }
        
        std::cerr << "Deff(pH=6)=" << get_alpha( pow(10.0,-6)) << std::endl;
        
        
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
