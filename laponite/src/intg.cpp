#include "yocto/math/fcn/functions.hpp"
#include "yocto/program.hpp"
#include "yocto/ios/ocstream.hpp"

using namespace yocto;
using namespace math;

#if 0
double I1(double a)
{
    const double h = 0.5*a;
    return 0.5*(exp(h)/pow(h,h)) * (exp( gamma_log(h) ) - gamma_i(h,h));
}
#endif

double logInvI1(double a)
{
    const double h = 0.5*a;
    return log(2.0) - h + h * log(h) - log( exp( gamma_log(h) ) - gamma_i(h,h)  );
}

YOCTO_PROGRAM_START()
{
    {
        ios::wcstream fp("i1.dat");
        for(double a=0.01;a<=8;a+=0.01)
        {
            const double lnY1 = logInvI1(a);
            fp("%15g %15g %15g %15g\n",log(a),lnY1,a,exp(lnY1));
        }
    }
}
YOCTO_PROGRAM_END()

