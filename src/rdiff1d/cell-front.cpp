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
