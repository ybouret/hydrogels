#include "yocto/sequence/vector.hpp"
#include "yocto/math/core/tao.hpp"
#include "yocto/program.hpp"
#include "yocto/physics/constants.hpp"
#include "yocto/math/core/tridiag.hpp"
#include "yocto/code/ipower.hpp"
#include "yocto/ios/ocstream.hpp"

using namespace yocto;
using namespace math;

#if 0
static const double Theta = Y_R * 298; //!< J/mol
static const double kH    = 29.41*1e2; //!< m^3.Pa/mol = J/mol
static const double rho   = Theta/kH;

static const double Patm = 1e5; //!< Pa

double getPressure( double t )
{
    return 300.0 * 100; //!< Pa
}


#endif

static const double Rini    = 0.1 * 1e-3; //!< initial size
static const double D       = 2e-9; //!< m^2/s

static const double Rend    = 1 * 1e-3; //!< final initial size


static size_t alpha = 1;


static void save_profile(const array<double> &r,
                         const array<double> &C)
{
    ios::wcstream fp("prof.dat");
    for(size_t i=1;i<=r.size();++i)
    {
        fp("%g %g\n", r[i], C[i]);
    }
    fp("\n");
}

static void save_q(const array<double> &r,
                   const array<double> &C)
{
    ios::acstream fp("q.dat");
    const size_t N = r.size();
    for(size_t i=1;i<=N;++i)
    {
        fp("%g %g\n", r[i]/r[1], (C[i]-C[1])/(C[N]-C[1]));
    }
    fp("\n");
}

#include "yocto/string/conv.hpp"

YOCTO_PROGRAM_START()
{

    double dotR = 1e-3/60;
    if(argc>1)
    {
        dotR = strconv::to<double>(argv[1],"dotR");
    }
    const size_t N = 1000;

    vector<double>  r(N);
    vector<double>  C(N);
    TriDiag<double> M(N);
    vector<double>  rhs(N);

    ios::ocstream::overwrite("prof.dat");
    ios::ocstream::overwrite("q.dat");
    ios::ocstream::overwrite("grad.dat");

    // initialize
    for(size_t i=1; i<=N; ++i )
    {
        r[i] = Rini + ((i-1)*Rend)/(N-1);
        C[i] = 1;
    }

    C[1] = 0.5;

    save_profile(r,C);

    const double dt  = 1e-3;
    const double Ddt = D * dt;

    double t = 0;

    for(size_t iter=1;iter<=20000;++iter)
    {
        t = iter * dt;
        const double R    = r[1];
        const double dR   = dt * dotR;

        // contact
        M.b[1] = 1;
        rhs[1] = 0.5;

        // core
        for(size_t i=2;i<N;++i)
        {
            rhs[i] =  C[i];
            // first order coefficient
            const double tdr  = r[i+1]-r[i-1];
            const double fac1 = (dR * ipower(R/r[i],alpha) - Ddt * alpha/r[i])/tdr;
            M.c[i] =  fac1;
            M.a[i] = -fac1;

            // second order term
            const double dr = tdr/2;
            const double aa = Ddt/(dr*dr);
            M.b[i]  = 1+aa+aa;;
            M.a[i] -= aa;
            M.c[i] -= aa;
        }


        // end
        M.a[N] = M.c[N] = 0;
        M.b[N] = 1;
        rhs[N] = 1.0; //!< final conc

        //std::cerr << "M=" << M << std::endl;


        // solve new concs
        M.solve(C,rhs);



        // move grid
        for(size_t i=1;i<=N;++i)
        {
            r[i] += dR * ipower(R/r[i],alpha);
        }

        // compute gradient
        const double drC = D*(C[2]-C[1])/(r[2]-r[1]);

        if( 0 == (iter%100) )
        {
            {
                ios::acstream fp("grad.dat");
                fp("%g %g\n", t, drC);
            }
            
            save_q(r,C);
        }
    }
    
    save_profile(r,C);
    
    
    
    
}
YOCTO_PROGRAM_END()

