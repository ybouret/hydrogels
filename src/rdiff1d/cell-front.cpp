#include "cell.hpp"

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
            pos = x0 + (x1-x0) * (v-y0)/(y1-y0);
            return true;
        }
    }
    return false;
}
