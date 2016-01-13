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


class Area
{
public:
    explicit Area() throw()
    {
    }

    virtual ~Area() throw()
    {
    }

    double Compute(const double t, const array<double> &a )
    {
        double t0 = a[1];
        double dt = t-t0;
        return a[2]+a[3]*dt+a[4]*dt*dt;
    }

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Area);
};

YOCTO_PROGRAM_START()
{

    //const double fac = numeric<double>::two_pi * D;

    for(int argi=1;argi<argc;++argi)
    {

        //______________________________________________________________________
        //
        // Load Data
        //______________________________________________________________________
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
        std::cerr << "#data=" << n << std::endl;


        //______________________________________________________________________
        //
        // Fitting log(p) then getting dp/p
        //______________________________________________________________________

        vector<double> lnp(n);
        vector<double> lnp_f(n);
        vector<double> dot_lnp(n);
        for(size_t i=1;i<=n;++i)
        {
            lnp[i] = log( pres[i] );
        }


        GLS<double>::Samples    samples(1);
        GLS<double>::Sample    &lpn_sample = samples.append(tmx,lnp,lnp_f);
        GLS<double>::Function   poly = _GLS::Create<double,_GLS::Polynomial>();
        GLS<double>::Proxy      poly_px(poly,min_of<size_t>(n-1,3));

        array<double>          &pcoef = poly_px.a;
        vector<bool>            pused(pcoef.size(),true);
        vector<double>          pcerr(pcoef.size());

        samples.prepare(pcoef.size());
        _GLS::Polynomial<double>::Start(lpn_sample, pcoef);
        if(!samples.fit_with(poly, pcoef, pused, pcerr))
        {
            throw exception("couldn't polynomial fit log(pressure) in %s", filename.c_str());
        }

        std::cerr << "Poly:" << std::endl;
        GLS<double>::display(std::cerr,pcoef, pcerr);
        numeric<double>::function poly_fn(&poly_px, & GLS<double>::Proxy::Compute);

        for(size_t i=1;i<=n;++i)
        {
            dot_lnp[i] = samples.diff(poly_fn,tmx[i]);
        }

        {
            ios::wcstream fp("smp.dat");
            for(size_t i=1;i<=n;++i)
            {
                fp("%g %g %g %g\n", tmx[i], lnp[i], lnp_f[i], dot_lnp[i]);
            }
        }


        //______________________________________________________________________
        //
        // Fitting area
        //______________________________________________________________________
        Area AA;
        GLS<double>::Proxy    AreaPx( &AA, & Area::Compute, 4);

        array<double> &acoef = AreaPx.a;
        vector<bool>   aused( acoef.size(), true );
        vector<double> acerr( acoef.size() );

        double &slope = acoef[3];
        //double &start = acoef[1];
        slope = (area[n]-area[1])/(tmx[n]-tmx[1]);

        vector<double> area_f(n);
        samples.release();
        (void) samples.append(tmx,area,area_f);
        samples.prepare(acoef.size());


        //aused[3] = false;

        if(!samples.fit_with(AreaPx.F, acoef, aused, acerr))
        {
            throw exception("couldn't fit area in %s", filename.c_str());
        }
        std::cerr << "Area:" << std::endl;
        GLS<double>::display(std::cerr,acoef,acerr);

        {
            ios::wcstream fp("sma.dat");
            for(size_t i=1;i<=n;++i)
            {
                fp("%g %g %g\n", tmx[i], area[i], area_f[i]);
            }
        }
    }
}
YOCTO_PROGRAM_END()
