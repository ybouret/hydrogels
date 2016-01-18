#include "yocto/sequence/vector.hpp"
#include "yocto/math/core/tao.hpp"
#include "yocto/program.hpp"
#include "yocto/physics/constants.hpp"
#include "yocto/math/core/tridiag.hpp"
#include "yocto/code/ipower.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/round.hpp"
#include "yocto/math/fit/glsf.hpp"
#include "yocto/eta.hpp"
#include "yocto/duration.hpp"

using namespace yocto;
using namespace math;

static const double Theta = Y_R * 298; //!< J/mol
static const double kH    = 29.41*1e2; //!< m^3.Pa/mol = J/mol





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

static void save_profile(const string        &p_name,
                         const array<double> &r,
                         const array<double> &C)
{
    ios::wcstream fp(p_name);
    for(size_t i=1;i<=r.size();++i)
    {
        fp("%g %g\n", r[i], C[i]);
    }
    fp("\n");
}

static void save_q(const string        &q_name,
                   const array<double> &s,
                   const array<double> &Q,
                   const array<double> &Qf)
{
    ios::wcstream fp(q_name);
    const size_t N = s.size();
    for(size_t i=1;i<=N;++i)
    {
        fp("%g %g %g\n", s[i], Q[i], Qf[i]);
    }
    fp("\n");
}

#include "yocto/string/conv.hpp"
#include "yocto/lua/lua-state.hpp"
#include "yocto/lua/lua-config.hpp"
#include "yocto/lua/lua-maths.hpp"

YOCTO_PROGRAM_START()
{
    eta ETA;

    Lua::State VM;
    lua_State *L = VM();

    double rho   = Theta/kH;
    std::cerr << "rho0=" << rho << std::endl;

    // number of points
    if(argc>1)
    {
        Lua::Config::DoFile(L,argv[1]);
    }

    const size_t N = size_t(Lua::Config::Get<lua_Number>(L, "N"));
    if(N<3)
    {
        throw exception("invalid #points!");
    }

    const unsigned alpha = unsigned(Lua::Config::Get<lua_Number>(L, "alpha"));
    if(alpha!=1&&alpha!=2)
    {
        throw exception("invalid alpha=%u", alpha);
    }

    // contact value
    Lua::Function<double>     _Cstar(L,"Cstar");
    numeric<double>::function Cstar(_Cstar);

    const double D      = Lua::Config::Get<lua_Number>(L,"D");
    const double R_ini  = Lua::Config::Get<lua_Number>(L,"R_ini");
    const double R_end  = Lua::Config::Get<lua_Number>(L,"R_end");
    double dt           = Lua::Config::Get<lua_Number>(L,"dt");
    double dt_save      = Lua::Config::Get<lua_Number>(L,"dt_sav");
    const  double t_run = Lua::Config::Get<lua_Number>(L,"t_run");


    const double delta_r =( R_end-R_ini)/(N-1);
    const double dt_max  = Square(delta_r)/D;
    std::cerr << "desired_dt=" << dt     << std::endl;
    std::cerr << "maximum_dt=" << dt_max << std::endl;


    dt = min_of<double>(dt,dt_max);
    const size_t every  = simulation_save_every(dt, dt_save);
    const size_t nIter  = simulation_iter(t_run, dt, every);
    
    
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
    used[2] = true;
    QShape shape;
    GLS<double>::Function F(&shape, &QShape::Compute );



    // initialize
    for(size_t i=1; i<=N; ++i )
    {
        r[i] = R_ini + (i-1)*delta_r;
        C[i] = 1;
    }
    r[N] = R_end;

    C[1] = Cstar(0);

    const string suffix    = vformat("%u_C%.2f.dat",alpha,C[1]);
    const string p_name = "p" + suffix;
    const string q_name = "q" + suffix;
    const string j_name = "j" + suffix;
    const string f_name = "f" + suffix;
    ios::ocstream::overwrite(p_name);
    ios::ocstream::overwrite(q_name);
    ios::ocstream::overwrite(j_name);
    ios::ocstream::overwrite(f_name);


    save_profile(p_name,r,C);


    const double Ddt    = D * dt;

    double t = 0;
    std::cerr << "#iter=" << nIter << std::endl;
    ETA.reset();
    std::cerr.flush();
    for(size_t iter=1;iter<=nIter;++iter)
    {
        t = iter * dt;
        const double R    = r[1];

        // compute flux
        const double J    = D*(C[2]-C[1])/(r[2]-r[1]);

        // compute dRdt
        const double dotR = rho * J / C[1];
        const double dR   = dt * dotR;

        //______________________________________________________________________
        //
        // contact
        //______________________________________________________________________

        M.b[1] = 1;
        rhs[1] = Cstar(t); // C*/Cinf

        //______________________________________________________________________
        //
        // core
        //______________________________________________________________________
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

        //______________________________________________________________________
        //
        // end
        //______________________________________________________________________
        M.a[N] = M.c[N] = 0;
        M.b[N] = 1;
        rhs[N] = 1.0; //!< final conc

        //std::cerr << "M=" << M << std::endl;

        //______________________________________________________________________
        //
        // solve new concs
        //______________________________________________________________________
        M.solve(C,rhs);

        //______________________________________________________________________
        //
        // move grid
        //______________________________________________________________________
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
                ios::acstream fp(j_name);
                fp("%g %g %g\n", t, J, (Q[2]-Q[1])/(s[2]-s[1]) );
            }

            if(!samples.fit_with(F,aorg,used,aerr) )
            {
                throw exception("Couldn't fit Q");
            }
            //std::cerr << "aorg=" << aorg << std::endl;
            {
                ios::acstream fp(f_name);
                fp("%g", t);
                for(size_t i=1;i<=nvar;++i)
                {
                    fp(" %g", aorg[i]);
                }
                fp("\n");
            }
            save_q(q_name,s,Q,Qf);
            save_profile(p_name,r,C);
            ETA(iter/double(nIter));
            const duration _left(ETA.time_left);
            fprintf(stderr,"t=%8.3lf | eta:%02u:%02u:%05.2f    \r",t,_left.h,_left.m,_left.s);
            fflush(stderr);
        }
    }
    std::cerr << std::endl;
    
    
}
YOCTO_PROGRAM_END()

