#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/fs/local-fs.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/exception.hpp"
#include "yocto/math/sig/smoother.hpp"
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
        
        size_t           count = 0;
        smoother<double> smooth;
        extender<double> xtd(extend_odd);

        smooth.upper_range = sm_dt/2;
        smooth.lower_range = sm_dt/2;
        smooth.degree      = sm_dg;
        
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
            // Fit Area
            //__________________________________________________________________
            vector<double> Anum(3,0);
            vector<double> Aden(1,0);
            vector<double> Afit(N,0);
            {
                vector<bool> usedP( Anum.size(), true );
                vector<bool> usedQ( Aden.size(), true );
                
                usedP[1] = false;
                Anum[1]  = A[1];
                
                least_squares<double>::sample SampleA(t,A,Afit);
                SampleA.Pade(Anum,usedP,Aden,usedQ);
            }
            
            //__________________________________________________________________
            //
            // Build Auxiliary quantities
            //__________________________________________________________________
            vector<double> PA(N,0);
            vector<double> smPA(N,0);
            vector<double> dPAdt(N,0);
            vector<double> smA(N,0);
            vector<double> dAdt(N,0);
            vector<double> PA32(N,0);
            vector<double> smPA32(N,0);
            vector<double> dPA32dt(N,0);
            
            smooth(smA,t,A,xtd,dAdt);

            for(size_t i=1;i<=N;++i)
            {
                PA[i]   = Pfit[i] * smA[i];
                PA32[i] = Pfit[i] * pow(smA[i],1.5);
            }
            
            smooth(smPA,  t,PA,  xtd,dPAdt);
            smooth(smPA32,t,PA32,xtd,dPA32dt);
            
            {
                const double omega2 = 0.375;
                const double omega3 = 0.288;
                ios::ocstream fp(outname,false);
                for(size_t i=1;i<=N;++i)
                {
                    const double growth2d = dPAdt[i]/pow(dAdt[i],omega2);
                    const double growth3d = (dPA32dt[i]/sqrt(smA[i]))/pow(dAdt[i],omega3);
                    //                                    1     2     3     4        5        6        7      8         9        10
                    fp("%g %g %g %g %g %g %g %g %g %g\n", t[i], A[i], P[i], Afit[i], Pfit[i], dAdt[i], PA[i], growth2d, PA32[i], growth3d);
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