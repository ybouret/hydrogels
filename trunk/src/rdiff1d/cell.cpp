#include "cell.hpp"

Cell:: ~Cell() throw()
{
}

Cell:: Cell( lua_State *L ) :
Collection(L),
ChemSys(L,*this),
Params(L,*this),
Workspace(vlayout,F,noghost),
M( Collection::size() ),
ini_side(L,"ini_side", *this),
ini_core(L,"ini_core", *this),
on_side(*this),
in_core(*this),
X( mesh.X() ),
pH( static_cast<Workspace&>(*this)["pH"].as<Array1D>() ),
pConc(M,as_capacity),
pFlux(M,as_capacity),
pIncr(M,as_capacity)
{
    //__________________________________________________________________________
    //
    // extra initialization
    //__________________________________________________________________________
    ini_side(*this,*this,0);
    on_side.load(C);
    std::cerr << "on_side=" << on_side << std::endl;
    
    ini_core(*this,*this,0);
    in_core.load(C);
    std::cerr << "in_core=" << in_core << std::endl;
    
    for(size_t i=0; i <= volumes;++i)
    {
        mesh.X()[i] = (i * length)/volumes;
    }
    
    Workspace &W = *this;
    for( Params::iterator i= Params::begin(); i != Params::end(); ++i )
    {
        const SpeciesParam &sp = **i;
        {
            Array1D &a = W[ sp.name ].as<Array1D>();
            pConc.push_back( &a );
        }
        
        {
            Array1D &a = W[ sp.flux ].as<Array1D>();
            pFlux.push_back( &a );
        }
        
        {
            Array1D &a = W[ sp.incr ].as<Array1D>();
            pIncr.push_back( &a );
        }
    }
    
    for( collection::iterator i = Collection::begin(); i != Collection::end(); ++i )
    {
        const species &sp = **i;
        const double   sD = sp.data.as<SpeciesData>().D;
        std::cerr << "D_" << sp.name << " = " << sD << std::endl;
        D.push_back( sD );
    }
    
}


void Cell:: init_all() throw()
{
    //std::cerr << "-- init all" << std::endl;
    for(size_t k=M;k>0;--k)
    {
        Array1D &c = *pConc[k];
        c[0] = on_side[k];
        const double value = in_core[k];
        for(size_t i=1;i<=volumes;++i)
            c[i] = value;
    }
    norm_all(0.0);
}


void Cell:: norm_all( double t )
{
    //std::cerr << "-- norm all" << std::endl;
    // assuming side is/sides are OK
    for(size_t i=1;i<=volumes;++i)
    {
        // fill in ChemSys.C
        for(size_t k=M;k>0;--k)
        {
            Array1D &c = *pConc[k];
            C[k] = c[i];
        }
        
        // chemistry
        normalize_C(t);
        
        // fill in workspace
        for(size_t k=M;k>0;--k)
        {
            Array1D &c = *pConc[k];
            c[i] = C[k];
        }
    }
    
    const Workspace &self = *this;
    const Array1D   &h    = self["H+"].as<Array1D>();
    for(size_t i=0;i<=volumes;++i)
    {
        pH[i] = -log10(h[i]);
    }
    
}

void Cell:: compute_fluxes()
{
    //std::cerr << "-- fluxes" << std::endl;
    
    for(size_t k=M;k>0;--k)
    {
        const Array1D &c  = *pConc[k];
        Array1D       &f  = *pFlux[k];
        const double   Dk = D[k];
        for(size_t i=0;i<volumes;++i)
        {
            f[i] = -Dk*(c[i+1]-c[i])/(X[i+1]-X[i]);
        }
        f[volumes]  = 0;
    }
    
}

#include "yocto/code/utils.hpp"

double Cell:: compute_increases( double dt, double t)
{
    //std::cerr << "-- increases " << std::endl;
    double shrink = -1;
    
    // core increases only
    const double fac = 2*dt;
    for(size_t i=1;i<volumes;++i)
    {
        
        for(size_t k=M;k>0;--k)
        {
            const Array1D &f  = *pFlux[k];
            Array1D       &I  = *pIncr[k];
            I[i]  = - fac * ( f[i]-f[i-1] )/(X[i+1]-X[i-1]);
            dC[k] = I[i];
        }
        
        legalize_dC(t);
        
        for(size_t k=M;k>0;--k)
        {
            Array1D       &I  = *pIncr[k];
            const Array1D &c  =  *pConc[k];
            const double   CC = c[i];
            const double   DD = I[i] = dC[k];
            if( DD<0 && -DD > CC)
            {
                const double slower = CC/(-DD);
                if(shrink<0)
                {
                    shrink = slower;
                }
                else
                {
                    shrink = min_of(shrink,slower);
                }
            }
        }
    }
    return shrink;
}

#include "yocto/ios/ocstream.hpp"

void Cell:: save_xy(const string &filename) const
{
    ios::ocstream fp(filename,false);
    Params::const_iterator j = Params::begin();
    for(size_t k=1;k<=M;++k,++j)
    {
        const SpeciesParam &sp = **j;
        const Array1D      &c  =  *pConc[k];
        fp("#[%s]\n", sp.name.c_str());
        for(size_t i=0;i<=volumes;++i)
        {
            fp("%g %g\n", X[i], c[i]);
        }
        fp("\n");
        
#if 0
        const Array1D &f = *pFlux[k];
        
        fp("#{%s}\n", sp.flux.c_str());
        for(size_t i=0;i<=volumes;++i)
        {
            fp("%g %g\n", X[i], f[i]);
        }
        fp("\n");
        
        const Array1D &I = *pIncr[k];
        fp("#<%s>\n", sp.incr.c_str());
        for(size_t i=0;i<=volumes;++i)
        {
            fp("%g %g\n", X[i], I[i]);
        }
        fp("\n");
#endif
    }
    
    fp("#pH\n");
    for(size_t i=0;i<=volumes;++i)
    {
        fp("%g %g\n", X[i], pH[i]);
    }

    
}


void Cell:: update_all(double factor)
{
    for(size_t i=1;i<volumes;++i)
    {
        
        for(size_t k=M;k>0;--k)
        {
            const Array1D &I  = *pIncr[k];
            Array1D       &c  = *pConc[k];
            c[i] += factor * I[i];
        }
        
    }
    
    // side
    for(size_t k=M;k>0;--k)
    {
        Array1D       &c  = *pConc[k];
        c[volumes] = c[volumes-1];
    }

    
}

void Cell::  step(double dt, double t)
{
    // we start from a valid conc
    const double t_end = t+dt;
    for(;;)
    {
        compute_fluxes();
        double shrink = compute_increases(dt, t);
        if(shrink>0)
        {
            const double fac = shrink/2;
            update_all(fac);
            dt *= fac;
            t += dt;
            dt = t_end - t;
            norm_all(t);
        }
        else
        {
            if(shrink<0)
            {
                update_all(1.0);
                norm_all(t_end);
                return;
            }
            else
            {
                throw exception("Invalid increases !");
            }
        }
    }
    
}

double Cell:: dx_min() const throw()
{
    double ans = fabs(X[1] - X[0]);
    
    for(size_t i=2;i<=volumes;++i)
    {
        const double tmp = fabs(X[i]-X[i-1]);
        if(tmp<ans) ans = tmp;
    }
    
    return ans;
}

double Cell:: D_max() const throw()
{
    double ans = 0;
    for(size_t i=1;i<=D.size(); ++i )
    {
        if(D[i]>ans) ans = D[i];
    }
    return ans;
}

