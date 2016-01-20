#include "yocto/program.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/math/sig/smooth.hpp"
#include "yocto/math/fit/glsf-spec.hpp"
#include "yocto/math/fcn/composition.hpp"
#include "yocto/math/point2d.hpp"
#include "yocto/sort/unique.hpp"

using namespace yocto;
using namespace math;

static  double D      = 2.2e-9; //!< m^2/s

const double phi2   = 1.525;
const double omega2 = 0.375;

const double phi3   = 2.308;
const double omega3 = 0.288;

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
            den += Square(a[4]) * t;
        }
        return num/den;
    }

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Area);
};

typedef point2d<double> coord_t;
static inline int compare_coord( const coord_t &lhs, const coord_t &rhs )
{
    return __compare<double>(lhs.x, rhs.x);
}


static double tscale = 0;
static double tstart = 0;


YOCTO_PROGRAM_START()
{

    GLS<double>::Samples samples(1);

    //__________________________________________________________________________
    //
    //
    // First Pass:: collecting all the pressure for the same experiment
    //
    //__________________________________________________________________________
    vector<coord_t> coords;


    for(int argi=1;argi<argc;++argi)
    {

        //______________________________________________________________________
        //
        // Load Data
        //______________________________________________________________________
        const string filename = argv[argi];
        vector<double> tt, pp;
        {
            ios::icstream fp(filename);
            data_set<double> ds;
            ds.use(1, tt);  // seconds
            ds.use(2, pp);  // mbar
            ds.load(fp);
        }
        const size_t n = tt.size();
        std::cerr << "#data=" << n << std::endl;
        for(size_t i=1;i<=n;++i)
        {
            coord_t c(tt[i],pp[i]*100);
            coords.push_back(c);
        }
    }

    std::cerr << "unique of " << coords.size() << std::endl;
    quicksort(coords,compare_coord);
    unique(coords,compare_coord);
    const size_t np = coords.size();
    std::cerr << "np=" << np<< std::endl;

    if(np<=3)
    {
        throw exception("Not Enough points...");
    }

    //__________________________________________________________________________
    //
    //
    // Building the d lnp / dt function
    //
    //__________________________________________________________________________
    vector<double> all_time(np);
    vector<double> all_pres(np);
    for(size_t i=1;i<=np;++i)
    {
        all_time[i] = coords[i].x;
        all_pres[i] = coords[i].y;
    }




    const double   pscale = tao::RMS(all_pres);
    tstart = all_time[1];
    tscale = all_time[np]-tstart;
    vector<double> red_all_time(np);
    vector<double> red_all_ln_p(np);
    vector<double> fit_all_ln_p(np);

    for(size_t i=1;i<=np;++i)
    {
        red_all_time[i] = (all_time[i]-tstart)/tscale;
        red_all_ln_p[i] = log( all_pres[i] / pscale );
    }

    {
        ios::wcstream fp("allpress.dat");
        for(size_t i=1;i<=np;++i)
        {
            fp("%g %g %g %g\n",all_time[i],all_pres[i],red_all_time[i],red_all_ln_p[i]);
        }
    }

    GLS<double>::Function     Poly = _GLS::Create<double,_GLS::Polynomial>();
    GLS<double>::Proxy        PolyPx(Poly,min_of<size_t>(np-1,4));
    numeric<double>::function PolyFn( &PolyPx, &GLS<double>::Proxy::Compute);
    array<double>            &pcoef = PolyPx.a;
    vector<bool>              pused(pcoef.size(),true);
    vector<double>            pcerr(pcoef.size());
    if( ! samples.fit_with(PolyPx.F, red_all_time, red_all_ln_p,  fit_all_ln_p, pcoef, pused, pcerr) )
    {
        throw exception("cannot fit log(p)");
    }
    std::cerr << "log(pressure):" << std::endl;
    GLS<double>::display(std::cerr, pcoef, pcerr);


    {
        ios::wcstream fp("fitpress.dat");
        for(size_t i=1;i<=np;++i)
        {
            fp("%g %g %g %g\n",all_time[i],red_all_ln_p[i],fit_all_ln_p[i], samples.diff(PolyFn,red_all_time[i])/tscale );
        }

    }

    {
        ios::wcstream fp("approx.dat");
        for(size_t i=2;i<np;++i)
        {
            const double dt   = all_time[i+1]-all_time[i-1];
            const double dlnp = log(all_pres[i+1]) - log(all_pres[i-1]);
            fp("%g %g\n",all_time[i],dlnp/dt);
        }

    }


    //__________________________________________________________________________
    //
    //
    // First Pass:: Building approx for areas
    //
    //__________________________________________________________________________
    ios::ocstream::overwrite("dota.dat");
    for(int argi=1;argi<argc;++argi)
    {

        //______________________________________________________________________
        //
        // Load Data
        //______________________________________________________________________
        const string filename = argv[argi];
        const string rootname = vfs::get_base_name(filename);
        vector<double> tmx, area, pres;
        {
            std::cerr << "Loading " << filename << std::endl;
            ios::icstream fp(filename);
            data_set<double> ds;
            ds.use(1, tmx);   // seconds
            ds.use(2, pres);  // mbar
            ds.use(3, area);  // m^2
            ds.load(fp);
        }
        const size_t n = tmx.size();
        std::cerr << "n=" << n << std::endl;
        if(n<=3)
        {
            throw exception("not enough in %s!", filename.c_str());
        }

        //______________________________________________________________________
        //
        // reduced area for fit
        //______________________________________________________________________
        vector<double> red_time(n);
        vector<double> red_area(n);
        vector<double> fit_area(n);
        const double ascale = tao::RMS(area);
        std::cerr << "ascale=" << ascale << std::endl;
        for(size_t i=1;i<=n;++i)
        {
            red_time[i] = (tmx[i]-tmx[1])/tscale;
            red_area[i] = area[i]/ascale;
            pres[i]    *= 100.0;
        }


        Area AA;
        GLS<double>::Proxy    AreaPx( &AA, & Area::Compute, min_of<size_t>(n-1,4) );

        array<double> &acoef = AreaPx.a;
        vector<bool>   aused( acoef.size(), true );
        vector<double> acerr( acoef.size() );

        samples.release();
        GLS<double>::Sample &AS = samples.append(red_time,red_area,fit_area);
        samples.prepare(acoef.size());

        {
            vector<double> atmp(3);
            _GLS::Polynomial<double>::Start(AS,atmp);
            for(size_t i=1;i<=3;++i) acoef[i] = atmp[i];
        }

        if(acoef.size()>=4)
        {
            acoef[4] = 0.01;
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
            fit_area[i] *= ascale;
            dot_area[i]  = ascale*samples.diff(AreaFn,red_time[i])/tscale;
        }

        {
            string outname = rootname;
            vfs::change_extension(outname,"fit.dat");
            std::cerr << "saving in " << outname << std::endl;
            ios::wcstream fp(outname);
            ios::acstream fp2("dota.dat");
            for(size_t i=1;i<=n;++i)
            {
                fp("%g %g %g %g\n", tmx[i], area[i], fit_area[i], dot_area[i]);
                fp2("%g %g\n", tmx[i], dot_area[i]);
            }
            fp2("\n");
        }

#if 0
        {
            const size_t NX = 1000;
            string outname = rootname;
            vfs::change_extension(outname,"xp.dat");
            ios::wcstream fp(outname);
            for(size_t i=0;i<=NX;++i)
            {
                const double t    = all_time[1] + (double(i)*(all_time[np]-all_time[1]))/NX;
                const double tred = (t - tmx[1])/tscale;
                fp("%g %g\n",t,ascale*AreaFn(tred));
            }
        }
#endif

        //______________________________________________________________________
        //
        //
        // Fitting curves
        //
        //______________________________________________________________________

        vector<double> y2(n);
        vector<double> y3(n);
        vector<double> dlnp(n);
        vector<double> invp(n);
        const double fac = numeric<double>::two_pi * D;
        for(size_t i=1;i<=n;++i)
        {
            const double tred = (tmx[i] - tstart)/tscale;
            invp[i] = 1.0/pres[i];
            dlnp[i] = samples.diff(PolyFn,tred)/tscale;
            const double tau     = area[i]     / fac;
            const double dot_tau = dot_area[i] / fac;

            y2[i]   = (dot_tau +             tau * dlnp[i]);
            y2[i]  /= phi2*pow(dot_tau,omega2);

            y3[i]   = (dot_tau + (2.0/3.0) * tau * dlnp[i]);
            y3[i]  /= phi3*pow(dot_tau,omega3);

        }

        {
            string outname = rootname;
            vfs::change_extension(outname,"curv.dat");
            std::cerr << "saving in " << outname << std::endl;
            ios::wcstream fp(outname);
            for(size_t i=1;i<=n;++i)
            {
                fp("%g %g %g\n",invp[i],y2[i],y3[i]);
            }
        }


    }



#if 0
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
    
    
}
#endif

}
YOCTO_PROGRAM_END()
