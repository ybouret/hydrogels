#include "yocto/program.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/math/sig/smooth.hpp"


using namespace yocto;
using namespace math;

YOCTO_PROGRAM_START()
{
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

        vector<double> sm_pres(n);
        vector<double> sm_area(n);
        vector<double> sm_pres_diff(n);
        vector<double> sm_area_diff(n);

        sm(xp,sm_area,tmx,area,sm_area_diff);
        sm(xp,sm_pres,tmx,pres,sm_pres_diff);

        {
            ios::wcstream fp("sm.dat");
            for(size_t i=1;i<=n;++i)
            {
                fp("%g %g %g %g %g\n",tmx[i],sm_pres[i],sm_area[i],sm_pres_diff[i],sm_area_diff[i]);
            }
        }



    }
}
YOCTO_PROGRAM_END()
