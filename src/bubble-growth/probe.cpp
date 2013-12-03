#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/fs/vfs.hpp"
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
        
        const double sm_dt  = strconv::to<double>(argv[1],"dt");
        const size_t sm_dg  = strconv::to<size_t>(argv[2],"degree");
        
        for( int k=3; k < argc; ++k )
        {
            //______________________________________________________________________
            //
            // Load raw data
            //______________________________________________________________________
            vector<double> t; //s
            vector<double> A; //mm2
            vector<double> P; //mbar
            {
                data_set<double> ds;
                ds.use(3, t);
                ds.use(4, A);
                ds.use(7, P);
                ios::icstream fp( argv[k] );
                ds.load(fp);
            }
            std::cerr << "Loaded " << t.size()  << " points" << std::endl;
            
            //______________________________________________________________________
            //
            // Remove bad data
            //______________________________________________________________________
            
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
            
            //______________________________________________________________________
            //
            // Build Auxiliary Qtty
            //______________________________________________________________________
            const size_t N = t.size();
            vector<double> AP(N,0);
            
            for(size_t i=1; i <= N; ++i )
            {
                A[i] *= 1.0e-6;      // SI
                P[i] *= 1.0e2;       // SI
                AP[i] = A[i] * P[i]; // SI
            }
            
            //______________________________________________________________________
            //
            // Build Smooth Quantity
            //______________________________________________________________________
            extend<double> xtd(extend_odd);
            
#if 0
            vector<double> rms_t;
            vector<double> rms_v;
            
            
            xtd.build_rms(rms_t, rms_v, t, A, degree);
            {
                ios::ocstream fp("rmsA.dat",false);
                for(size_t i=1; i <= rms_t.size();++i)
                {
                    fp("%g %g\n", rms_t[i], rms_v[i]);
                }
            }
            
            xtd.build_rms(rms_t, rms_v, t, P, degree);
            {
                ios::ocstream fp("rmsP.dat",false);
                for(size_t i=1; i <= rms_t.size();++i)
                {
                    fp("%g %g\n", rms_t[i], rms_v[i]);
                }
            }
            
            xtd.build_rms(rms_t, rms_v, t, AP, degree);
            {
                ios::ocstream fp("rmsAP.dat",false);
                for(size_t i=1; i <= rms_t.size();++i)
                {
                    fp("%g %g\n", rms_t[i], rms_v[i]);
                }
            }
#endif
            
            
            vector<double> smA(N,0);
            vector<double> smdAdt(N,0);
            vector<double> smP(N,0);
            xtd(smA,t,A,sm_dt,sm_dg,&smdAdt);
            xtd(smP,t,P,sm_dt,sm_dg,0);
            
            
            //______________________________________________________________________
            //
            // Build Auxiliary Qtty
            //______________________________________________________________________
            
            std::cerr << "Loaded " << N  << " points" << std::endl;
            {
                ios::ocstream fp("bubble.dat",false);
                for(size_t i=1;i<=N;++i)
                {
                    fp("%g %g %g %g %g %g %g\n", t[i], A[i], P[i], AP[i], smA[i], smdAdt[i], smP[i]);
                }
            }
            
            
        }
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