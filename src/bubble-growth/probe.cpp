#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/fs/local-fs.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/exception.hpp"
#include "yocto/math/sig/smoother.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/math/fit/least-squares.hpp"
#include "yocto/math/polynomial.hpp"
#include "yocto/code/utils.hpp"

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
        
        ios::ocstream::overwrite("slopes.dat");
        
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
            // Fit Pressure with a line
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
            // Fit Area with Pade (3,1), fixing A[1]
            //__________________________________________________________________
            
            vector<double> Afit(N,0);
            vector<double> dAfitdt(N,0);
            
            {
                vector<double> Anum(3,0);
                vector<double> Aden(1,0);
                vector<bool> usedP( Anum.size(), true );
                vector<bool> usedQ( Aden.size(), true );
                
                usedP[1] = false;
                Anum[1]  = A[1];   // ??
                
                least_squares<double>::sample SampleA(t,A,Afit);
                SampleA.Pade(Anum,usedP,Aden,usedQ);
                
                polynomial<double> PP;
                for(size_t i=1;i<=Anum.size();++i) PP.add(i-1,Anum[i]);
                polynomial<double> QQ; QQ.add(0,1);
                for(size_t i=1;i<=Aden.size();++i) QQ.add(i,Aden[i]);
                
                std::cerr << "R(x)=" << PP << "/" << QQ << std::endl;
                
                //______________________________________________________________
                //
                // Deduce dAfitdt
                //______________________________________________________________
                polynomial<double>::derivative(PP,QQ);
                std::cerr << "R'(x)=" << PP << "/" << QQ << std::endl;
                for(size_t i=1;i<=N;++i)
                {
                    dAfitdt[i] = PP(t[i])/QQ(t[i]);
                }
            }
            
            //__________________________________________________________________
            //
            // Compute growth factors
            //__________________________________________________________________
            vector<double> PA(N,0);
            vector<double> PA32(N,0);
            for(size_t i=1;i<=N;++i)
            {
                PA[i]   = Pfit[i] * Afit[i];
                PA32[i] = Pfit[i] * pow(Afit[i],1.5);
            }
            
            vector<double> dPAdt(N,0);
            vector<double> dPA32dt(N,0);
            
            for(size_t i=2;i<N;++i)
            {
                dPAdt[i]   = (PA[i+1]  -PA[i-1]  )/(t[i+1]-t[i-1]);
                dPA32dt[i] = (PA32[i+1]-PA32[i-1])/(t[i+1]-t[i-1]);
            }
            dPAdt[1] = (4*PA[2]   - 3*PA[1] - PA[3]  )/(t[3]  -t[1]);
            dPAdt[N] = (4*PA[N-1] - 3*PA[N] - PA[N-2])/(t[N-2]-t[N]);
            
            
            dPA32dt[1] = (4*PA32[2]   - 3*PA32[1] - PA32[3]  )/(t[3]  -t[1]);
            dPA32dt[N] = (4*PA32[N-1] - 3*PA32[N] - PA32[N-2])/(t[N-2]-t[N]);
            
            vector<double> g2d(N,0);
            vector<double> g3d(N,0);
            
            //const double tpd    = 2 * M_PI * (2e-9);
            const double omega2 = 0.375;
            const double omega3 = 0.288;
            const double Phi2   = 1.525;
            const double Phi3   = 2.308;
            
            const double Theta  = 8.3144621 * (273+20);
            
            for(size_t i=1;i<=N;++i)
            {
               
                const double rate = dAfitdt[i];
                const double rho2 = Theta * pow(2*M_PI,(1-omega2))*Phi2*pow(rate,omega2);
                const double rho3 = Theta * 3.0/pow(2*M_PI,omega3)*Phi3*pow(rate,omega3);
                
                g2d[i] = dPAdt[i];
                g3d[i] = dPA32dt[i]/sqrt(Afit[i]);
                
                g3d[i] /= rho2;
                g3d[i] /= rho3;
            }

        
            
            
            ios::ocstream fp(outname,false);
            for(size_t i=1;i<=N;++i)
            {
                //                                         1    2    3    4       5       6          7     8        9       10         11     12
                fp("%g %g %g %g %g %g %g %g %g %g %g %g\n",t[i],A[i],P[i],Afit[i],Pfit[i],dAfitdt[i],PA[i],dPAdt[i],PA32[i],dPA32dt[i],g2d[i],g3d[i]);
            }
            
            
            const size_t   ns = min_of<size_t>(16,N);
            vector<double> xx(ns,0);
            vector<double> yy(ns,0);
            vector<double> zz(ns,0);
            
            least_squares<double>::sample fit_slope(xx,yy,zz);
            
            const size_t   dof = 2;
            vector<double> a2d(dof,0);
            vector<double> e2d(dof,0);
            vector<bool>   uad(dof,true);
            for(size_t i=1;i<=ns;++i)
            {
                xx[i] = Pfit[i];
                yy[i] = g2d[i];
                zz[i] = 0;
            }
            fit_slope.polynomial(a2d, uad, e2d);
            
            vector<double> a3d(dof,0);
            vector<double> e3d(dof,0);
            for(size_t i=1;i<=ns;++i)
            {
                xx[i] = Pfit[i];
                yy[i] = g3d[i];
                zz[i] = 0;
            }
            fit_slope.polynomial(a3d, uad, e3d);
            
            
            const double D = 2e-9;
            const double D2d = pow(D,1-omega2);
            const double D3d = pow(D,1-omega3);
            const double slope2d = a2d[2];
            const double slope3d = a3d[2];
            {
                ios::ocstream fp("slopes.dat",true);
                fp("%u %g %g %g %g #%s\n",unsigned(count+1),slope2d,slope3d,-D2d/slope2d,-D3d/slope3d,vfs::get_base_name(argv[k]));
            }
            
            
            
#if 0
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
#endif
            
            
            
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