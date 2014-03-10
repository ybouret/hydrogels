#include "yocto/math/types.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/exception.hpp"
#include "yocto/ios/rc.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ptr/auto.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/math/dat/spline.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/fcn/zfind.hpp"
#include "yocto/string/conv.hpp"


using namespace yocto;
using namespace math;;


class DGel : public spline1D<double>
{
public:
    explicit DGel( const vector<double> &C, const vector<double> &Rho ) :
    spline1D<double>(spline_natural,C,Rho),
    match(0.5)
    {
    }
    
    double match;
    
    virtual ~DGel() throw() {}
    
    double ZeroFunc( double c )
    {
        return get(c) - match;
    }
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(DGel);
};

// %gel
// 1          2          4          8
// ratios
// 0.322881   0.226271   0.144068   0.076271
int main(int argc, char *argv[] )
{
    try
    {
        
        // prepare the chemical slow down
        vector<double> C;
        vector<double> Rho;
        
        {
            ios::resources         rc( argv[0] );
            auto_ptr<ios::istream> fp( rc.load_stream("Dgel.out"));
            data_set<double> ds;
            ds.use(1,C);
            ds.use(3,Rho);
            ds.load(*fp);
        }
        
        DGel gel(C,Rho);
        
#if 0
        const size_t N = C.size();
        ios::ocstream fp("dpline.dat",false);
        for(size_t i=0; i<=10000;++i)
        {
            const double cc = C[1] + (i*(C[N]-C[1])/10000.0);
            const double dd = gel(cc);
            fp("%g %g\n", cc, dd );
        }
#endif
        
        if(argc>1)
        {
            vector<double> percent;
            vector<double> ratios;
            
            {
                data_set<double> ds;
                ds.use(1,percent);
                ds.use(2,ratios);
                ios::icstream fp(argv[1]);
                ds.load(fp);
            }
            
            const size_t nr = ratios.size();
            std::cerr << "percent= " << percent << std::endl;
            std::cerr << "ratios = " << ratios  << std::endl;
            numeric<double>::function zfunc( &gel, & DGel::ZeroFunc);
            
            zfind<double> solve( 1e-8 );
            const size_t N = 100;
            for(size_t r=1;r<=nr;++r)
            {
                const double per = percent[r];
                const double rho = ratios[r];
                if(rho<=0)
                    throw exception("invalid ratio %g", rho );
                const double sqrt_rho = sqrt(rho);
                const double phi_min  = (sqrt_rho+sqrt_rho)/(1.0+sqrt_rho);
                std::cerr << "rho    = " << rho << std::endl;
                std::cerr << "phi_min= " << phi_min << std::endl;
                const double phi_max  = 1;
                
                
                
                const string fn = vformat("rho_gel%.5gp.dat",per);
                std::cerr << "--> <" << fn <<  ">" << std::endl;
                
                gel.match = rho;
                const double Cmax     = solve(zfunc,0,C.back());
                std::cerr << "Cmax= " << Cmax << std::endl;
                
                ios::ocstream fp(fn,false);
                fp("#phi%g Cg\n",rho);
                for( size_t i=0; i<=N; ++i )
                {
                    const double phi = i <= 0 ? phi_min : (i>=N ? phi_max : phi_min + (i * (phi_max-phi_min)/N ) );
                    const double sqrt_fac = phi/(2.0-phi);
                    const double fac      = sqrt_fac * sqrt_fac;
                    gel.match             = rho/fac;
                    double Cg = 0;
                    if(i>0)
                    {
                        Cg = solve(zfunc,0,C.back());
                    }
                    //std::cerr << phi << " " << gel.match << " " << Cg << std::endl;
                    fp("%.5e %.5e\n", phi, Cg );
                }
                
            }
            
#if 0
            vector<double> poro;
            poro.push_back(1);
            poro.push_back(0.95);
            poro.push_back(0.90);
            poro.push_back(0.85);
            poro.push_back(0.80);
            poro.push_back(0.75);
            poro.push_back(0.70);
            
            std::cerr << " ratio";
            for(size_t j=1; j<=poro.size(); ++j)
            {
                std::cerr << " " << poro[j];
            }
            std::cerr << std::endl;
            for(size_t i=1; i <= ratios.size(); ++i )
            {
                const double ratio = ratios[i];
                std::cerr << ratio;
                for(size_t j=1; j <= poro.size();++j)
                {
                    const double phi      = poro[j];
                    const double sqrt_fac = phi/(2.0-phi);
                    const double fac      = sqrt_fac * sqrt_fac;
                    gel.match             = ratio/fac;
                    const double Cgel     = solve(zfunc,0,C.back());
                    std::cerr << " " << -log10(Cgel);
                }
                std::cerr << std::endl;
            }
#endif
            
        }
        
        
        
        
        return 0;
    }
    catch(const exception &e )
    {
        std::cerr << e.what() << std::endl;
        std::cerr << e.when() << std::endl;
    }
    catch(...)
    {
    }
    return -1;
}