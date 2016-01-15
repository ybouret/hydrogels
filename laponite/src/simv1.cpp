#include "yocto/sequence/vector.hpp"
#include "yocto/math/core/tao.hpp"
#include "yocto/program.hpp"
#include "yocto/physics/constants.hpp"
#include "yocto/math/core/tridiag.hpp"
#include "yocto/code/ipower.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/round.hpp"
#include "yocto/math/fit/glsf.hpp"

using namespace yocto;
using namespace math;

static const double Theta = Y_R * 298; //!< J/mol
static const double kH    = 29.41*1e2; //!< m^3.Pa/mol = J/mol



static const double Rini    = 0.1 * 1e-3; //!< initial size
static const double D       = 2e-9; //!< m^2/s

static const double Rend    = 1 * 1e-3; //!< final initial size


static size_t alpha = 1;


class QShape
{
public:
    inline QShape() {}
    inline ~QShape() throw() {}

    double Compute(const double s, const array<double> &a )
    {
        const double X = s-1;
        double       P = 0;
        for(size_t i=1;i<=a.size();++i)
        {
            P += a[i] * ipower(X,i);
        }
        return 1.0 - exp( -P );
    }


private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(QShape);
};

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

static void save_q(const array<double> &s,
                   const array<double> &Q,
                   const array<double> &Qf)
{
    ios::wcstream fp("q.dat");
    const size_t N = s.size();
    for(size_t i=1;i<=N;++i)
    {
        fp("%g %g %g\n", s[i], Q[i], Qf[i]);
    }
    fp("\n");
}

#include "yocto/string/conv.hpp"

YOCTO_PROGRAM_START()
{

    double rho   = Theta/kH;
    std::cerr << "rho0=" << rho << std::endl;

    if(argc>1)
    {
        alpha = strconv::to<size_t>(argv[1],"alpha");
        if(alpha!=1&&alpha!=2)
        {
            throw exception("invalid alpha!");
        }
    }

    const size_t N = 1000;

    vector<double>  r(N);
    vector<double>  C(N);
    TriDiag<double> M(N);
    vector<double>  rhs(N);
    vector<double>  s(N);
    vector<double>  Q(N);
    vector<double>  Qf(N);


    GLS<double>::Samples samples(1);
    samples.append(s, Q, Qf);

    const size_t nvar = 3;
    vector<double> aorg(nvar);
    vector<bool>   used(nvar,true);
    vector<double> aerr(nvar);
    samples.prepare(nvar);
    aorg[1] = 1;
    used[2] = false;
    QShape shape;
    GLS<double>::Function F(&shape, &QShape::Compute );


    ios::ocstream::overwrite("prof.dat");
    ios::ocstream::overwrite("q.dat");
    ios::ocstream::overwrite("grad.dat");
    ios::ocstream::overwrite("coef.dat");

    // initialize
    for(size_t i=1; i<=N; ++i )
    {
        r[i] = Rini + ((i-1)*Rend)/(N-1);
        C[i] = 1;
    }

    C[1] = 0.5;

    save_profile(r,C);

    double dt           = 1e-4;
    double dt_save      = 0.1;
    const  double t_run = 20;
    const size_t every  = simulation_save_every(dt, dt_save);
    const size_t nIter  = simulation_iter(t_run, dt, every);
    const double Ddt    = D * dt;

    double t = 0;
    std::cerr << "#iter=" << nIter << std::endl;
    for(size_t iter=1;iter<=nIter;++iter)
    {
        t = iter * dt;
        const double R    = r[1];

        // compute flux
        // compute gradient
        const double J    = D*(C[2]-C[1])/(r[2]-r[1]);
        const double dotR = rho * J / C[1];
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


        if( 0 == (iter%every) )
        {

            for(size_t i=1;i<=N;++i)
            {
                s[i] = r[i]/r[1];
                Q[i] = (C[i]-C[1])/(C[N]-C[1]);
            }

            {
                ios::acstream fp("grad.dat");
                fp("%g %g %g\n", t, J, (Q[2]-Q[1])/(s[2]-s[1]) );
            }

            if(!samples.fit_with(F,aorg,used,aerr) )
            {
                throw exception("Couldn't fit Q");
            }
            //std::cerr << "aorg=" << aorg << std::endl;
            {
                ios::acstream fp("coef.dat");
                fp("%g", t);
                for(size_t i=1;i<=nvar;++i)
                {
                    fp(" %g", aorg[i]);
                }
                fp("\n");
            }
            save_q(s,Q,Qf);
            save_profile(r,C);
            std::cerr << "."; std::cerr.flush();
        }
    }
    std::cerr << std::endl;
    

}
YOCTO_PROGRAM_END()

