#include "yocto/sequence/vector.hpp"
#include "yocto/math/types.hpp"
#include "yocto/math/fcn/intg.hpp"
#include "yocto/ios/ocstream.hpp"

using namespace yocto;
using namespace math;

namespace
{
    class wrapper
    {
    public:
        explicit wrapper() :
        lam(1),
        call_f2(this, & wrapper::compute_f2),
        call_f3(this, & wrapper::compute_f3)
        {
        }
        
        virtual ~wrapper() throw()
        {
        }
        
        double             lam;
        integrator<double> intg;
        
        inline double get_f2(double Lambda)
        {
            lam = Lambda;
            double ans = 0;
            for(size_t i=1;i<200;++i)
            {
                ans += intg(i,i+1,call_f2,1e-7);
            }
            return ans;
        }
        
        inline double get_f3(double Lambda)
        {
            lam = Lambda;
            double ans = 0;
            for(size_t i=1;i<200;++i)
            {
                ans += intg(i,i+1,call_f3,1e-7);
            }
            return ans;
        }
        

        
    private:
        YOCTO_DISABLE_COPY_AND_ASSIGN(wrapper);
        
        inline double compute_f2(double s)
        {
            return pow(s,lam-1.0) * exp(-lam*0.5*(s*s-1.0));
        }
        
        inline double compute_f3(double s)
        {
            return exp(lam*(1.5 - 1.0/s -0.5 *s*s))/(s*s);
        }
        
        numeric<double>::function call_f2;
        numeric<double>::function call_f3;

        
        
    };
}

int main(int argc, char *argv[] )
{
    
    try
    {
        const int    top = 50;
        const double rho = 10;
        wrapper w;
        const double alpha = log10(rho);
        
        
        const double I21 = w.get_f2(1);
        std::cerr << "I21=" << I21 << std::endl;
        
        if(false)
        {
            ios::ocstream fp("I2.dat",false);
            for(int i=-top;i<=top;++i)
            {
                const double lam = pow(10.0,(i*alpha)/top);
                std::cerr << "lam=" << lam << std::endl;
                const double I2 = w.get_f2(lam);
                fp("%g %g %g %g\n",lam,I2,log(lam),log(I2));
            }
        }
        
        const double I31 = w.get_f3(1);
        std::cerr << "I31=" << I31 << std::endl;
        if(true)
        {
            ios::ocstream fp("I3.dat",false);
            for(int i=-top;i<=top;++i)
            {
                const double lam = pow(10.0,(i*alpha)/top);
                std::cerr << "lam=" << lam << std::endl;
                const double I3 = w.get_f3(lam);
                fp("%g %g %g %g\n",lam,I3,log(lam),log(I3));
            }
        }
        
    }
    catch(...)
    {
        std::cerr << "unhandled error" << std::endl;
    }
    return 1;
}