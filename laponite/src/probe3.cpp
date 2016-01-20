#include "yocto/program.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/math/sig/smooth.hpp"
#include "yocto/math/fit/glsf-spec.hpp"



using namespace yocto;
using namespace math;

static  double D      = 2.0e-9; //!< m^2/s


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

    ios::ocstream::overwrite("coef.dat");

    // for fit
    GLS<double>::Samples    samples(1);

    // for smoothing
    smooth<double>          sm;
    sm.degree      = 3;
    sm.upper_range = 0.5;
    sm.lower_range = 0.5;

    // for expansion
    expand<double>          xp(expand_odd);


    for(int argi=1;argi<argc;++argi)
    {

        //______________________________________________________________________
        //
        // Load Data
        //______________________________________________________________________
        const string filename = argv[argi];
        const string rootname = vfs::get_base_name(filename);
        std::cerr << "rootname=" << rootname << std::endl;
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
        // preparing reduced data for fitting
        //______________________________________________________________________
        vector<double> red_time(n);
        vector<double> red_area(n);
        vector<double> invp(n);
        vector<double> red_invp(n);

        for(size_t i=1;i<=n;++i)
        {
            pres[i] *= 100;             //-- pressure in pascal
            invp[i]  = 1.0 / pres[i];   //-- Pa^(-1)
        }

        const double tscale = tmx[n] - tmx[1];
        const double ascale = tao::RMS(area);
        const double iscale = tao::RMS(invp);

        std::cerr << "tscale=" << tscale << std::endl;
        std::cerr << "ascale=" << ascale << std::endl;
        std::cerr << "iscale=" << iscale << std::endl;

        for(size_t i=1;i<=n;++i)
        {
            red_area[i] = area[i] / ascale;
            red_time[i] = (tmx[i]-tmx[1]) / tscale;
            red_invp[i] = invp[i] / iscale;
        }

        //______________________________________________________________________
        //
        // fitting invp
        //______________________________________________________________________
        vector<double> fit_invp(n);


        GLS<double>::Function     Poly = _GLS::Create<double,_GLS::Polynomial>();
        GLS<double>::Proxy        PolyPx(Poly,min_of<size_t>(n-1,3));
        numeric<double>::function PolyFn( &PolyPx, &GLS<double>::Proxy::Compute);
        array<double>          &pcoef = PolyPx.a;
        vector<bool>            pused(pcoef.size(),true);
        vector<double>          pcerr(pcoef.size());
        if( ! samples.fit_with(PolyPx.F, red_time, red_invp, fit_invp, pcoef, pused, pcerr) )
        {
            throw exception("cannot fit invp for %s", filename.c_str() );
        }

        std::cerr << "invp:" << std::endl;
        GLS<double>::display(std::cerr, pcoef, pcerr);

        for(size_t i=1;i<=n;++i)
        {
            fit_invp[i] *= iscale;
        }

        {
            string outname = rootname;
            vfs::change_extension(outname, "pres.fit.dat");
            std::cerr << "outname=" << outname << std::endl;
            ios::wcstream fp(outname);
            for(size_t i=1;i<=n;++i)
            {
                fp("%g %g %g\n", tmx[i], invp[i], fit_invp[i]);
            }

        }


#if 0
        //______________________________________________________________________
        //
        // smooth...
        //______________________________________________________________________
        sm(xp, lnpr_f, tred, lnpr, dot_lnp);
        for(size_t i=1;i<=n;++i)
        {
            dot_lnp[i] = pscale * dot_lnp[i] / tscale;
        }

        {
            string outname = rootname;
            vfs::change_extension(outname, "pres.fit.dat");
            std::cerr << "outname=" << outname << std::endl;
            ios::wcstream fp(outname);
            for(size_t i=1;i<=n;++i)
            {
                fp("%g %g %g %g\n", tmx[i], log(pres[i]), lnpr_f[i]*pscale, dot_lnp[i]);
            }
        }
#endif


#if 0
        //______________________________________________________________________
        //
        // Fitting lnpr then deducing d(log(p))/dt
        //______________________________________________________________________
        GLS<double>::Function     Poly = _GLS::Create<double,_GLS::Polynomial>();
        GLS<double>::Proxy        PolyPx(Poly,min_of<size_t>(n-1,3));
        numeric<double>::function PolyFn( &PolyPx, &GLS<double>::Proxy::Compute);
        array<double>          &pcoef = PolyPx.a;
        vector<bool>            pused(pcoef.size(),true);
        vector<double>          pcerr(pcoef.size());
        if( ! samples.fit_with(PolyPx.F, tred, lnpr, lnpr_f, pcoef, pused, pcerr) )
        {
            throw exception("cannot with ln(pressure) for %s", filename.c_str() );
        }

        std::cerr << "reduced ln(pressure):" << std::endl;
        GLS<double>::display(std::cerr, pcoef, pcerr);

        for(size_t i=1;i<=n;++i)
        {
            dot_lnp[i] = pscale * samples.diff(PolyFn,tred[i])/tscale;
        }

        {
            string outname = rootname;
            vfs::change_extension(outname, "pres.fit.dat");
            std::cerr << "outname=" << outname << std::endl;
            ios::wcstream fp(outname);
            for(size_t i=1;i<=n;++i)
            {
                fp("%g %g %g %g\n", tmx[i], log(pres[i]), lnpr_f[i]*pscale, dot_lnp[i]);
            }
        }
#endif


#if 0
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
            pres[i] *= 100;            //-- pressure in pascal
            lnp[i]   = log( pres[i] ); //-- log(pressure/pascal)
            ipr[i]   = 1.0 / pres[i];  //-- Pa^(-1)
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
        // Fitting curves
        //______________________________________________________________________
        vector<double> tau(n);
        vector<double> dot_tau(n);
        vector<double> y2(n);
        vector<double> y3(n);
        const double fac = numeric<double>::two_pi * D;
        for(size_t i=1;i<=n;++i)
        {
            tau[i]     = area[i]/fac;
            dot_tau[i] = dot_area[i]/fac;
            y2[i]      = dot_tau[i]+ tau[i] * dot_lnp[i];
            y3[i]      = dot_tau[i]+ (2.0/3.0) * tau[i] * dot_lnp[i];
        }



        vector<double> y2f(n);
        vector<double> y3f(n);

        const size_t   nf = 2;
        vector<double> p2(nf), p2err(nf);
        vector<double> p3(nf), p3err(nf);
        vector<bool>   used(nf,true);

        if(!samples.fit_with(poly, ipr, y2, y2f, p2, used, p2err))
        {
            throw exception("2D fit failure");
        }
        else
        {
            std::cerr << "2D:" << std::endl;
            GLS<double>::display(std::cerr, p2, p2err);
        }



        if(!samples.fit_with(poly, ipr, y3, y3f, p3, used, p3err))
        {
            throw exception("3D fit failure");
        }
        else
        {
            std::cerr << "3D:" << std::endl;
            GLS<double>::display(std::cerr, p3, p3err);
        }

        
        
        
        
        {
            ios::wcstream fp("fit.dat");
            for(size_t i=1;i<=n;++i)
            {
                fp("%g %g %g %g %g\n", ipr[i], y2[i], y3[i], y2f[i], y3f[i]);
            }
        }
        
        {
            ios::acstream fp("coef.dat");
            fp("#%s\n", filename.c_str());
            fp("%d %g %g %g %g\n", argi, p2[2], p3[2], p2[1], p3[1]);
        }
        
#endif
        
    }
}
YOCTO_PROGRAM_END()
