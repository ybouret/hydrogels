#include "yocto/sequence/vector.hpp"
#include "yocto/math/types.hpp"
#include "yocto/math/fcn/intg.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/fit/glsf-spec.hpp"

using namespace yocto;
using namespace math;

namespace
{
    class wrapper
    {
    public:
        static const size_t IMAX = 200;

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
            for(size_t i=1;i<IMAX;++i)
            {
                ans += intg(i,i+1,call_f2,1e-7);
            }
            return 1.0/ans;
        }

        inline double get_f3(double Lambda)
        {
            lam = Lambda;
            double ans = 0;
            for(size_t i=1;i<IMAX;++i)
            {
                ans += intg(i,i+1,call_f3,1e-7);
            }
            return 1.0/ans;
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


#include "yocto/program.hpp"

YOCTO_PROGRAM_START()
{
    const double min_p10 = -1;
    const double max_p10 =  1;

    wrapper w;
    const double Y21 = w.get_f2(1);
    std::cerr << "Y21=" << Y21 << std::endl;

    const double Y31 = w.get_f3(1);
    std::cerr << "Y31=" << Y31 << std::endl;

    size_t N = 100;



    vector<double> lna(N);
    vector<double> ly2(N),ly2f(N);
    vector<double> ly3(N),ly3f(N);

    {
        ios::wcstream fp("y.dat");
        ios::wcstream fp2("ly.dat");
        for(size_t i=1;i<=N;++i)
        {
            const double p  = min_p10 + ( double(i-1) * (max_p10-min_p10) )/double(N-1);
            const double a  = pow(10.0,p);
            const double Y2 = w.get_f2(a);
            const double Y3 = w.get_f3(a);
            fp("%g %g %g\n", a, Y2, Y3);
            lna[i] = log(a);
            ly2[i] = log(Y2);
            ly3[i] = log(Y3);
            fp2("%g %g %g\n", lna[i], ly2[i], ly3[i] );
        }
    }

    GLS<double>::Function poly = _GLS::Create<double,_GLS::Polynomial>();

    size_t nvar=3;
    vector<double> p2(nvar),p2err(nvar),p3(nvar),p3err(nvar);

    p2[1] = log(Y21);
    p3[1] = log(Y31);
    vector<bool> used(nvar,true);
    used[1] = false;

    GLS<double>::Samples samples(1);

    if( !samples.fit_with(poly, lna, ly2, ly2f, p2, used, p2err) )
    {
        throw exception("Cannot fit 2D");
    }
    std::cerr << "2D:" << std::endl;
    GLS<double>::display(std::cerr, p2, p2err);

    if( !samples.fit_with(poly, lna, ly3, ly3f, p3, used, p3err) )
    {
        throw exception("Cannot fit 2D");
    }
    std::cerr << "3D:" << std::endl;
    GLS<double>::display(std::cerr, p3, p3err);

    {
        ios::wcstream fp("ly_fit.dat");
        for(size_t i=1;i<=N;++i)
        {
            fp("%g %g %g\n", lna[i], ly2f[i], ly2f[i]);
        }
    }
    
    
}
YOCTO_PROGRAM_END()