#include "cell.hpp"

double Cell:: max_dt() const
{
    double Dmax = 0;
    for( library::const_iterator i = begin(); i != end(); ++i )
    {
        const double D = (**i).get<SpeciesData>().D;
        if( D > Dmax ) Dmax = D;
    }
    if( Dmax <= 0 )
        throw exception("Invalid Diffusion Coefficients");
    return alpha * (dx*dx) / Dmax;
}


double Cell:: locate( const double pH ) const
{
    const WorkspaceBase &self = *this;
    const double h0 = pow(10.0,-pH);
    const Array &h  = self["H+"].as<Array>();
    for(unit_t i=X.lower;i<X.upper;++i)
    {
        const double hlo = h[i];
        const double hup = h[i+1];
        if( (h0-hlo)*(hup-h0) >= 0 )
        {
            // found
            const double pH_lo = -log10(hlo);
            const double pH_up = -log10(hup);
            const double num   = pH - pH_lo;
            const double den   = pH_up - pH_lo;
            if( pH_lo < pH_up )
            {
                if(pH<=pH_lo)
                    return X[i];
                else
                {
                    if(pH>=pH_up)
                    {
                        return X[i];
                    }
                    else
                        return X[i] + (dx * num) / den;
                }
            }
            else
            {
                if( pH_lo > pH_up )
                {
                    if(pH>=pH_lo)
                        return X[i];
                    else
                    {
                        if(pH<=pH_up)
                            return X[i+1];
                        else
                            return X[i] + (dx * num) / den;
                    }
                }
                else
                {
                    return 0.5 * (X[i]+X[i+1]);
                }
            }
            
        }
    }
    throw exception("Not found pH=%g", pH);
}