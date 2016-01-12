#include "yocto/program.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/math/sig/smooth.hpp"
#include "yocto/math/fit/glsf-spec.hpp"


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

    //const double fac = numeric<double>::two_pi * D;
    GLS<double>::Function poly = _GLS::Create<double,_GLS::Polynomial>();
    GLS<double>::Function pade = _GLS::Create<double,_GLS::Pade>();

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

        vector<double> ipr(n);
        vector<double> ipr_f(n);
        for(size_t i=1;i<=n;++i)
        {
            ipr[i] = 1.0/pres[i];
        }

        vector<double> pres_f(n);
        vector<double> area_f(n);

        //______________________________________________________________________
        //
        // Pressure
        //______________________________________________________________________
        GLS<double>::Samples samples(1);
        GLS<double>::Sample &ps = samples.append(tmx, ipr, ipr_f);

        vector<double> pa(3);
        samples.prepare(pa.size());
        _GLS::Polynomial<double>::Start(ps, pa);
        std::cerr << "predicted=" << pa << std::endl;

        {
            ios::wcstream fp("sm.dat");
            for(size_t i=1;i<=n;++i)
            {
                fp("%g %g %g\n", tmx[i], pres[i], 1.0/ipr_f[i] );
            }
        }


    }
}
YOCTO_PROGRAM_END()
