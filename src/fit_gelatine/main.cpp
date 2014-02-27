#include "yocto/math/types.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/exception.hpp"
#include "yocto/ios/rc.hpp"
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

int main(int argc, char *argv[] )
{
    try
    {
        
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
            numeric<double>::function zfunc( &gel, & DGel::ZeroFunc);
            
            zfind<double> solve( 1e-8 );
            
            vector<double> ratios;
            
            for(int i=1; i<argc;++i)
            {
                ratios.push_back( strconv::to<double>(argv[i], "ratio") );
            }
            
            vector<double> poro;
            poro.push_back(1);
            poro.push_back(0.95);
            poro.push_back(0.90);
            poro.push_back(0.85);
            poro.push_back(0.80);
            poro.push_back(0.75);
            //poro.push_back(0.70);
            
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