#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/fs/local-fs.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/exception.hpp"
#include "yocto/math/sig/extend.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/math/fit/least-squares.hpp"

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
            
            {
                const double t0 = t[1];
                for(size_t i=1; i <= N;++i)
                {
                    t[i] -= t0;
                }
            }
            
            
            for(size_t i=1; i <= N; ++i )
            {
                A[i] *= 1.0e-6;      // SI: m^2
                P[i] *= 1.0e2;       // SI: Pa
            }
            
            
            
            //__________________________________________________________________
            //
            // Fit Pressure
            //__________________________________________________________________
            vector<double> Pfit(N,0);
            vector<double> Pcof(2,0);
            vector<double> Perr(2,0);
            {
                vector<bool> used(2,true);
                least_squares<double>::sample SampleP(t,P,Pfit);
                SampleP.polynomial(Pcof, used, Perr);
            }
            std::cerr << "Pcof=" << Pcof << " +/- " << Perr << std::endl;
            const double Pdot = Pcof[2];
            std::cerr << "slope=" << Pdot << std::endl;
            
            //__________________________________________________________________
            //
            // Build Auxiliary quantities
            //__________________________________________________________________
            extend2<double> xtd(extend_odd);
            vector<double>  smA(N,0);
            vector<double>  dAdt(N,0);
            vector<double>  AP(N,0);
            vector<double>  smAP(N,0);
            vector<double>  dAPdt(N,0);

            xtd(smA,t,A,sm_dt,sm_dg,dAdt); // smooth A and compute dAdt
            for(size_t i=1;i<=N;++i)
            {
                AP[i]   = Pfit[i] * smA[i];
            }
            xtd(smAP,t,AP,sm_dt,sm_dg,dAPdt); // smooth AP -> dAPdt
            
            {
                ios::ocstream fp(outname,false);
#if 0
                fp("#t A P Pfit smA dAdt AP dAPdt dPdt mech\n");
                for(size_t i=1;i<=N;++i)
                {
                    const double mech = - Pdot * smA[i] / Pfit[i];
                    const double diff = dAPdt[i]/Pfit[i];
                    //                                       1     2     3     4        5       6        7       8        9     10          11
                    fp("%g %g %g %g %g %g %g %g %g %g %g\n", t[i], A[i], P[i], Pfit[i], smA[i], dAdt[i], AP[i], dAPdt[i], Pdot, mech, diff);
                }
#endif
                const double omega2 = 0.375;
                const double omega3 = 0.288;
                fp("#t A P smA Pfit dAdt dPdt PdAdt AdPdt growth2D growthd3D\n");
                for(size_t i=1;i<=N;++i)
                {
                    const double PdAdt = Pfit[i] * dAdt[i];
                    const double AdPdt = smA[i]  * Pdot;
                    const double g2d   = (PdAdt + AdPdt)/pow(dAdt[i],omega2);
                    const double g3d   = (PdAdt + 2.0/3.0 * AdPdt)/pow(dAdt[i],omega3);
                    
                    //                                       1     2     3     4        5        6        7     8      9      10   11
                    fp("%g %g %g %g %g %g %g %g %g %g %g\n", t[i], A[i], P[i], smA[i],  Pfit[i], dAdt[i], Pdot, PdAdt, AdPdt, g2d, g3d );
                }
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