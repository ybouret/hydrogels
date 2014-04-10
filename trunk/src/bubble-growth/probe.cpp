#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/fs/local-fs.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/exception.hpp"
#include "yocto/math/sig/extend.hpp"
#include "yocto/string/conv.hpp"

#include <iostream>

using namespace yocto;
using namespace math;



int main(int argc, char *argv[])
{
    
    const char *prog = vfs::get_base_name(argv[0]);
    try
    {
        if(argc<=3)
            throw exception("usage: %s dt degree [datafiles]",prog);
        
        const double sm_dt  = strconv::to<double>(argv[1],"dt");        // smoothing window size
        const size_t sm_dg  = strconv::to<size_t>(argv[2],"degree");    // smoothing polynomial degree
        
        vfs &fs = local_fs::instance();
        string  outdir = "pbubbles";
        fs.as_directory(outdir);
        std::cerr << "will save in [" << outdir << "]" << std::endl;
        fs.create_sub_dir(outdir);
        fs.remove_files_with_extension_in(outdir, "dat");
        
        size_t count = 0;
        for( int k=3; k < argc; ++k )
        {
            //__________________________________________________________________
            //
            // Load raw data
            //__________________________________________________________________
            vector<double> t; //s
            vector<double> A; //mm2
            vector<double> P; //mbar
            std::cerr << std::endl << "Loading [" << argv[k] << "]" << std::endl;
            {
                data_set<double> ds;
                ds.use(3, t);
                ds.use(4, A);
                ds.use(7, P);
                ios::icstream fp( argv[k] );
                ds.load(fp);
            }
            std::cerr << "Loaded " << t.size()  << " points" << std::endl;
            string outname = outdir + vfs::get_base_name(argv[k]);
            vfs::change_extension(outname, "dat");
            
            //__________________________________________________________________
            //
            // Remove bad data
            //__________________________________________________________________
            
            while( t.size() && t.front() <= 0 )
            {
                t.pop_front();
                A.pop_front();
                P.pop_front();
            }
            
            while(t.size() && t.back() <= 0 )
            {
                t.pop_back();
                A.pop_back();
                P.pop_back();
            }
            
            
            
            //__________________________________________________________________
            //
            // Build Auxiliary Qtty
            //__________________________________________________________________
            const size_t N = t.size();
            if(N<=0)
            {
                std::cerr << "No Points" << std::endl;
                continue;
            }
            std::cerr << "#Effective Points = " << N << std::endl;
            //std::cerr << "t=" << t << std::endl;
            
            {
                const double t0 = t[1];
                for(size_t i=1; i <= N;++i)
                {
                    t[i] -= t0;
                }
            }
            vector<double> AP(N,0);
            
            for(size_t i=1; i <= N; ++i )
            {
                A[i] *= 1.0e-6;      // SI: m^2
                P[i] *= 1.0e2;       // SI: Pa
                AP[i] = A[i] * P[i]; // SI...
            }
            
            //__________________________________________________________________
            //
            // Build Smooth Quantity
            //__________________________________________________________________
            extend<double> xtd(extend_odd);
            
            
            
            vector<double> smA(N,0);            // smoothed area
            vector<double> smdAdt(N,0);         // smoothed area time derivative
            vector<double> smP(N,0);            // smoothed pressure
            vector<double> smdPdt(N,0);         // smoothed pressure derivative
            vector<double> smAP(N,0);           // smoothed area x pressure
            vector<double> smdAPdt(N,0);        // smoothed (area x pressure) time derivative
            std::cerr << "\tsmoothing Area" << std::endl;
            xtd(smA,t,A,sm_dt,sm_dg,&smdAdt);
            std::cerr << "\tsmoothing Pressure" << std::endl;
            xtd(smP,t,P,sm_dt,sm_dg,&smdPdt);
            std::cerr << "\tsmoothing Area x Pressure" << std::endl;
            xtd(smAP,t,AP,sm_dt,sm_dg,&smdAPdt);
            
            //__________________________________________________________________
            //
            // Build Auxiliary Qtty
            //__________________________________________________________________
            std::cerr << "Loaded " << N  << " points" << std::endl;
            
#if 0
            {
                ios::ocstream fp("bubble.dat",false);
                fp << "#t A P AP smA smdAdt smP smdPdt smAP smdAPdt\n";
                for(size_t i=1;i<=N;++i)
                {
                    //                                    1     2     3     4      5       6          7       8          9        10
                    fp("%g %g %g %g %g %g %g %g %g %g\n", t[i], A[i], P[i], AP[i], smA[i], smdAdt[i], smP[i], smdPdt[i], smAP[i], smdAPdt[i] );
                }
            }
            
            {
                ios::ocstream fp("dt_area.dat",false);
                for(size_t i=1;i<=N;++i)
                {
                    const double mech    = - A[i] * smdPdt[i] / P[i];
                    const double sm_mech = - smA[i] * smdPdt[i] / smP[i];
                    const double diff    = smdAPdt[i] / P[i];
                    const double sm_diff = smdAPdt[i] / smP[i];
                    //                                   1     2     3              4          5                 6     7       8     9
                    fp("%g %g %g %g %g %g %g %g %g\n", t[i], A[i], smdAdt[i], (mech+diff), (sm_mech+sm_diff), mech, sm_mech, diff, sm_diff);
                }
                
            }
#endif
            
            //__________________________________________________________________
            //
            // Build Auxiliary Qtty
            //__________________________________________________________________
            {
                ios::ocstream fp(outname,false);
                
            }
            ++count;
        }
        std::cerr << std::endl;
        std::cerr << "processed " << count << " files" << std::endl;
        return 0;
    }
    catch( const exception &e )
    {
        std::cerr << "*** in " << prog     << std::endl;
        std::cerr << "*** "    << e.what() << std::endl;
        std::cerr << "*** "    << e.when() << std::endl;
    }
    catch(...)
    {
        std::cerr << "*** unhandled exception in " << prog << std::endl;
    }
    return -1;
}