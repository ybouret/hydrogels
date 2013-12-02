#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/fs/vfs.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/exception.hpp"
#include "yocto/math/sig/extend.hpp"

#include <iostream>

using namespace yocto;
using namespace math;


int main(int argc, char *argv[])
{
    
    const char *prog = vfs::get_base_name(argv[0]);
    try
    {
        if(argc<=1)
            throw exception("need a data file");
        
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
            ios::icstream fp( argv[1] );
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
        vector<double> F(N,0);
        vector<double> dF(N,0);
        
        for(size_t i=1; i <= N; ++i )
        {
            A[i] *= 1.0e-6; // SI
            P[i] *= 1.0e2;  // SI
            AP[i] = A[i] * P[i];
        }
        
        extend<double> xtd(extend_odd);
        
        
        //______________________________________________________________________
        //
        // Build Auxiliary Qtty
        //______________________________________________________________________
        
        std::cerr << "Loaded " << N  << " points" << std::endl;
        {
            ios::ocstream fp("bubble.dat",false);
            for(size_t i=1;i<=N;++i)
            {
                fp("%g %g %g %g %g %g\n", t[i], A[i], P[i], AP[i], F[i], dF[i]);
            }
        }
        return 0;
    }
    catch(...)
    {
        std::cerr << "unhandled exception in " << prog << std::endl;
    }
    return -1;
}