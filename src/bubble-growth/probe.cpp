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
            if(N<=0)
            {
                std::cerr << "No Points" << std::endl;
                continue;
            }
            
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
                A[i] *= 1.0e-6;      // SI
                P[i] *= 1.0e2;       // SI
                AP[i] = A[i] * P[i]; // SI
            }
            
            //______________________________________________________________________
            //
            // Build Smooth Quantity
            //______________________________________________________________________
            extend<double> xtd(extend_odd);
            
            
            
            vector<double> smA(N,0);
            vector<double> smdAdt(N,0);
            vector<double> smP(N,0);
            vector<double> smdPdt(N,0);
            vector<double> smAP(N,0);
            vector<double> smdAPdt(N,0);
            xtd(smA,t,A,sm_dt,sm_dg,&smdAdt);
            xtd(smP,t,P,sm_dt,sm_dg,&smdPdt);
            xtd(smAP,t,AP,sm_dt,sm_dg,&smdAPdt);
            
            //______________________________________________________________________
            //
            // Build Auxiliary Qtty
            //______________________________________________________________________
            
            std::cerr << "Loaded " << N  << " points" << std::endl;
            {
                ios::ocstream fp("bubble.dat",false);
                fp << "#t A P AP smA smdAdt smP smdPdt smAP smdAPdt\n";
                for(size_t i=1;i<=N;++i)
                {
                    //                                    1     2     3     4      5       6          7       8          9        10
                    fp("%g %g %g %g %g %g %g %g %g %g\n", t[i], A[i], P[i], AP[i], smA[i], smdAdt[i], smP[i], smdPdt[i], smAP[i], smdAPdt[i] );
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