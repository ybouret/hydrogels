#include "cell.hpp"
#include "yocto/math/fcn/zfind.hpp"


double Cell:: zfront(double x) throw()
{
    return front_ip.polint(x) - search_value;
}

bool Cell:: find_front( double &pos ) const
{
    const Workspace &self = *this;
    const Array1D   &A    = self[search_field].as<Array1D>();
    const double     v    = search_value;
    for(size_t i=0,j=1; i<volumes;++i,++j)
    {
        const double x0 = X[i];
        const double x1 = X[j];
        const double y0 = A[i];
        const double y1 = A[j];
        
        // do we cross the y0 -> y1 segment ?
        if( (v-y0)*(v-y1) < 0 )
        {
            front_ip.free();
            front_ip.insert(x0, y0);
            front_ip.insert(x1, y1);
            {
                if(i>0)
                {
                    const size_t im=i-1;
                    front_ip.insert(X[im],A[im]);
                }
                if(j<volumes)
                {
                    const size_t jp=j+1;
                    front_ip.insert(X[jp],A[jp]);
                }
                front_ip.optimize();
            }
            zfind<double> solve( 0 );
            pos = solve(front_fn,x0,x1);
            return true;
        }
    }
    return false;
}

static inline
double get_curvature( const Array1D &X, const Array1D &Y, unit_t i )
{
    const double x0 = X[i];
    const double f0 = Y[i];
    
    const double xm = X[i-1] - x0;
    const double fm = Y[i-1];
    
    const double xp = X[i+1] - x0;
    const double fp = Y[i+1];
    
    const double xm2 = xm*xm;
    const double xp2 = xp*xp;
    
    const double num = xp*(fm-f0) - xm*(fp-f0);
    
    //std::cerr << "f=" << fm << " " << f0 << " " << fp << std::endl;
    return num / (xm2*xp-xp2*xm);
    
}



#include "yocto/code/utils.hpp"

bool Cell:: find_inflection(double &xx, double &yy) const
{
    const Workspace &self = *this;
    const Array1D   &A    = self[search_field].as<Array1D>();
    double x0 = X[2];
    double y0 = A[2];
    double c0 = get_curvature(X,A,2);
    
    const unit_t itop = X.upper - 1;
    for(unit_t i=3;i<itop;++i)
    {
        const double x1 = X[i];
        const double y1 = A[i];
        const double c1 = get_curvature(X,A,i);
        const double dc = c1-c0;
        if( c0 * c1 <= 0.0 && ( fabs(dc)>0 ) )
        {
            //std::cerr << "c0=" << c0 << std::endl;
            //std::cerr << "c1=" << c1 << std::endl;
            
            const double dx = x1-x0;
            
            xx = clamp<double>(x0,x0 - c0*dx/dc,x1);
            yy = y0 + (xx-x0)*(y1-y0)/dx;
            
            return true;
        }
        
        x0 = x1;
        y0 = y1;
        c0 = c1;
    }
    
    return false;
}

