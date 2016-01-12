#include "yocto/program.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/math/sig/smooth.hpp"


using namespace yocto;
using namespace math;

static  double D      = 2.0e-8; //!< m^2/s
static  double omega2 = 0.375;
static  double phi2   = 1.525;
static  double omega3 = 0.288;
static  double phi3   = 2.208;


double F2D(double dotlam,
           double corlam)
{
    return (dotlam + corlam) / (phi2*pow(dotlam,omega2));
}

double F3D(double dotlam,
           double corlam)
{
    return (dotlam + (2.0/3.0) * corlam) / (phi3*pow(dotlam,omega3));
}


YOCTO_PROGRAM_START()
{

    const double fac = numeric<double>::two_pi * D;

    for(int argi=1;argi<argc;++argi)
    {
        const string filename = argv[argi];
        vector<double> tmx, pres, area;

        {
            ios::icstream fp(filename);
            data_set<double> ds;
            ds.use(1, tmx);
            ds.use(2, pres);
            ds.use(3, area);
            ds.load(fp);
        }
        const size_t n = tmx.size();
        const double dt_full = tmx[n] - tmx[1];
        std::cerr << "#data=" << n << std::endl;

        const expand<double> xp(expand_odd);
        smooth<double>       sm;
        const double dt = dt_full/2;
        sm.degree      = 5;
        sm.lower_range = dt/2;
        sm.upper_range = dt/2;

        vector<double> smx_pres(n);
        vector<double> smx_area(n);
        vector<double> dot_pres(n);
        vector<double> dot_area(n);
        vector<double> lnp(n);
        vector<double> smx_lnp(n);
        vector<double> dot_lnp(n);
        for(size_t i=1;i<=n;++i)
        {
            lnp[i] = log(pres[i]);
        }

        sm(xp,smx_area,tmx,area,dot_area);
        sm(xp,smx_pres,tmx,pres,dot_pres);
        sm(xp,smx_lnp, tmx,lnp, dot_lnp);

        vector<double> lam(n);
        vector<double> dlm(n);
        vector<double> y2d(n);
        vector<double> y3d(n);
        vector<double> cor(n);
        vector<double> ipr(n);

        for(size_t i=1;i<=n;++i)
        {
            const double lambda = area[i]/fac;
            const double dotlam = dot_area[i]/fac;
            const double corlam = dot_lnp[i] * lambda;
            lam[i]  = lambda;
            dlm[i]  = dotlam;
            cor[i]  = corlam;
            y2d[i]  = F2D(dotlam, corlam);
            y3d[i]  = F3D(dotlam, corlam);
            ipr[i]  = 1.0/pres[i];
        }

        {
            ios::wcstream fp("sm.dat");
            for(size_t i=1;i<=n;++i)
            {
                fp("%g %g %g %g %g %g %g\n", tmx[i], ipr[i], y2d[i], y3d[i], lam[i], dlm[i], -cor[i]);
                //fp("%g %g %g %g %g %g %g %g %g %g\n", tmx[i], area[i], smx_area[i], dot_area[i]/fac, pres[i], smx_pres[i], dot_pres[i]/pres[i], lnp[i], smx_lnp[i], dot_lnp[i]);
            }
        }



    }
}
YOCTO_PROGRAM_END()
