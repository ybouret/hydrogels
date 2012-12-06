#include "yocto/math/fcn/zfind.hpp"
#include "yocto/math/fcn/functions.hpp"

#include "yocto/string/vfs-utils.hpp"
#include "yocto/exception.hpp"
#include "yocto/ios/ocstream.hpp"

using namespace yocto;

static const double Kw  =  1e-14;
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

static inline double d(double t, double x)
{
    if(x<=0)
        return d0;
    else
    {
		const double arg = x/sqrt(4*D*t);
        return d1 + (d0-d1) * math::qerfc( x/sqrt(4*D*t) );
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

static inline double h(double t, double x )
{
    return 0.5 * ( d(t,x) + s(t,x) );
}

int main(int argc, char *argv[])
{
    const char *progname = _vfs::get_base_name(argv[0]);
    try
    {
        
        double t  = 5;
        double dx = 0.00001;
        {
			ios::ocstream fp( "profil.dat", false );
	        for( double x=dx; x <= 0.002; x += dx )
	        {
				fp("%g %g\n", x, -log10(h(t,x)) );
			}
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
