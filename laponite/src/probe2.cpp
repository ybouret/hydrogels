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
        double num = a[1] + a[2] * t + a[3] * t*t;
        double den = 1;
        if(a.size()>=4)
        {
            den += a[4] * t;
        }
        return num/den;
    }

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Area);
};

YOCTO_PROGRAM_START()
{


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
            ds.use(1, tmx);   // seconds
            ds.use(2, pres);  // mbar
            ds.use(3, area);  // m^2
            ds.load(fp);
        }
        const size_t n = tmx.size();
        std::cerr << "#data=" << n << std::endl;
        if(n<=3)
        {
            throw exception("not enough points!");
        }

        //______________________________________________________________________
        //
        // Fitting log(p) then getting dp/p
        //______________________________________________________________________

        vector<double> lnp(n);
        vector<double> lnp_f(n);
        vector<double> dot_lnp(n);
        vector<double> ipr(n);
        for(size_t i=1;i<=n;++i)
        {
            pres[i] *= 100;
            lnp[i]   = log( pres[i] );
            ipr[i]   = 1.0 / pres[i];
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

        vector<double> ared(n); //!< reduced area for units
        vector<double> afit(n); //!< fitted value
        vector<double> tred(n); //!< reduced time

        double ascale = tao::RMS(area);
        double tscale = tmx[n] - tmx[1];
        std::cerr << "ascale=" << ascale << std::endl;
        for(size_t i=1;i<=n;++i)
        {
            ared[i] = area[i]/ascale;
            tred[i] = (tmx[i]-tmx[1])/tscale;
        }



        Area AA;
        GLS<double>::Proxy    AreaPx( &AA, & Area::Compute, min_of<size_t>(n-1,4) );

        array<double> &acoef = AreaPx.a;
        vector<bool>   aused( acoef.size(), true );
        vector<double> acerr( acoef.size() );

        samples.release();
        GLS<double>::Sample &AS = samples.append(tred,ared,afit);
        samples.prepare(acoef.size());

        {
            vector<double> atmp(3);
            _GLS::Polynomial<double>::Start(AS,atmp);
            for(size_t i=1;i<=3;++i) acoef[i] = atmp[i];
        }

        if(!samples.fit_with(AreaPx.F, acoef, aused, acerr))
        {
            throw exception("couldn't fit area in %s", filename.c_str());
        }
        std::cerr << "Area:" << std::endl;
        GLS<double>::display(std::cerr,acoef,acerr);

        numeric<double>::function AreaFn( &AreaPx, & GLS<double>::Proxy::Compute);

        vector<double> dot_area(n);
        for(size_t i=1;i<=n;++i)
        {
            afit[i]    *= ascale;
            dot_area[i] = ascale*samples.diff(AreaFn,tred[i])/tscale;
        }

        {
            ios::wcstream fp("sma.dat");
            for(size_t i=1;i<=n;++i)
            {
                fp("%g %g %g %g\n", tmx[i], area[i], afit[i], dot_area[i]);
            }
        }

        //______________________________________________________________________
        //
        // Constructing curves
        //______________________________________________________________________
        const double fac = numeric<double>::two_pi * D;
        vector<double> lambda(n);
        vector<double> dotlam(n);
        vector<double> corlam(n);
        vector<double> y2(n);
        vector<double> y3(n);
        for(size_t i=1;i<=n;++i)
        {
            lambda[i] = area[i]/fac;
            dotlam[i] = dot_area[i]/fac;
            corlam[i] = dot_lnp[i] * lambda[i];
            y2[i]     = F2D(dotlam[i],corlam[i]);
            y3[i]     = F3D(dotlam[i],corlam[i]);
        }

        {
            ios::wcstream fp("lam.dat");
            for(size_t i=1;i<=n;++i)
            {
                fp("%g %g %g %g\n", tmx[i],y2[i], y3[i], ipr[i]);
            }
        }

        //______________________________________________________________________
        //
        // Fitting...
        //______________________________________________________________________
        vector<double> p2(2), p3( p2.size() ), p2err( p2.size() ), p3err( p2.size() );
        vector<double> y2f(n);
        vector<double> y3f(n);

        vector<bool> lused(p2.size(),true);

        samples.release();
        _GLS::Polynomial<double>::Start( samples.append(ipr,y2,y2f), p2);
        samples.prepare(p2.size());
        if( ! samples.fit_with(poly,p2,lused,p2err) )
        {
            throw exception("unexpected fit failure for 2D");
        }

        std::cerr << "Fit2D: " << std::endl;
        GLS<double>::display(std::cerr, p2, p2err);


        samples.release();
        _GLS::Polynomial<double>::Start( samples.append(ipr,y3,y3f), p3);
        samples.prepare(p3.size());
        if( ! samples.fit_with(poly,p3,lused,p3err) )
        {
            throw exception("unexpected fit failure for 3D");
        }

        std::cerr << "Fit3D: " << std::endl;
        GLS<double>::display(std::cerr, p3, p3err);

        {
            ios::wcstream fp("fit.dat");
            for(size_t i=1;i<=n;++i)
            {
                fp("%g %g %g %g %g\n",ipr[i],y2[i],y2f[i],y3[i],y3f[i]);
            }
        }






    }
}
YOCTO_PROGRAM_END()
