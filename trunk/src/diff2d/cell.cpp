#include "cell.hpp"



Cell:: Cell( lua_State *L ) :
Library(L),
ChemSys(*this,L),
Parameters(*this,L),
Workspace( *this, *this ),
ini_bath( "ini_bath", *this, L),
ini_skin( "ini_skin", *this, L),
shrink(1)
{
    
    //--------------------------------------------------------------------------
    //-- initialize bath
    //--------------------------------------------------------------------------
    solution S( *this );
    ini_bath( *this, 0.0 );
    S.get( C );
    std::cerr << "*** bath=" << std::endl << S << std::endl;
    Coord u(0,0);
    for( ; u.x <= volumes.x; ++u.x )
    {
        saveC(C,u);
    }
    std::cerr << "*** \tbath is loaded..." << std::endl;
    
    //--------------------------------------------------------------------------
    //-- initialize
    //--------------------------------------------------------------------------
    ini_skin( *this, 0.0 );
    S.get( C );
    std::cerr << "*** skin=" << std::endl << S << std::endl;
    for( unit_t j=1; j <= volumes.y; ++j )
    {
        for( unit_t i=0; i <= volumes.x; ++i )
        {
            const Coord u(i,j);
            saveC(C,u);
        }
    }
    std::cerr << "*** \tskin is loaded..." << std::endl;
    
    porosity.ld(1);
    
    vtk.save("init.vtk", "init", *this, var, vlayout);
}


Cell:: ~Cell() throw()
{
    
}

#if 0
void Cell:: compute_fluxes()
{
    //size_t q = 0;
    for( Specs::iterator p = specs.begin(); p != specs.end(); ++p )
    {
        //++q;
        SpeciesData        &sd = **p; assert( sd.U ); assert( sd.F );
        const Array        &U  = *(sd.U);
        VertexArray        &F  = *(sd.F);
        for( unit_t j=vtop.y;j>0;--j)
        {
            F[j][0].ldz();
            F[j][volumes.x].ldz();
            for( unit_t i=vtop.x;i>0;--i)
            {
                const double D = sd.D * porosity[j][i];
                Vertex &g = F[j][i];
                g.x = D * (U[j][i+1]-U[j][i-1])/twodel.x;
                g.y = D * (U[j+1][i]-U[j-1][i])/twodel.y;
            }
        }
        F[0][0].ldz();         F[0][volumes.x].ldz();
        F[volumes.y][0].ldz(); F[volumes.y][volumes.x].ldz();
    }
}
#endif

void Cell:: compute_increases(double t, double dt)
{
#if 0
    for( unit_t j=vtop.y;j>0;--j)
    {
        for( unit_t i=vtop.x;i>0;--i)
        {
            size_t q = 0;
            for( Specs::iterator p = specs.begin(); p != specs.end(); ++p )
            {
                ++q;
                SpeciesData        &sd = **p;
                const Array        &U  = *(sd.U);
                const VertexArray  &F  = *(sd.F);
                Array              &I  = *(sd.I);
                C[q] = U[j][i];
                const double dvg = (F[j][i+1].x - F[j][i-1].x)/twodel.x + (F[j+1][i].y - F[j-1][i].y)/twodel.y;
                dC[q] = dvg*dt;
                
                //-- chemistry
                reduce(t);
                I[j][i] = dC[q];
            }
        }
    }
#endif
    for( unit_t j=vtop.y;j>0;--j)
    {
        const unit_t jm = j-1;
        const unit_t jp = j+1;
        
        for( unit_t i=vtop.x;i>0;--i)
        {
            const unit_t im = i-1;
            const unit_t ip = i+1;
            size_t q = 0;
            for( Specs::iterator p = specs.begin(); p != specs.end(); ++p )
            {
                ++q;
                SpeciesData        &sd = **p;
                const Array        &U  = *(sd.U);
                //const VertexArray  &F  = *(sd.F);
                Array              &I  = *(sd.I);
                const double U0 = C[q] = U[j][i];
                
                const double D0 = sd.D * porosity[j][i];
                
                const double Dyp = 0.5 * ( D0 + sd.D * porosity[jp][i]);
                const double Dym = 0.5 * ( D0 + sd.D * porosity[jm][i]);
                const double Iy  = inv_sq.y * (Dyp * ( U[jp][i] - U0 ) - Dym * (U0 - U[jm][i]) );
                
                const double Dxp = 0.5 * ( D0 + sd.D * porosity[j][ip]);
                const double Dxm = 0.5 * ( D0 + sd.D * porosity[j][im]);
                const double Ix  = inv_sq.x * (Dxp * ( U[j][ip] - U0 ) - Dxm * (U0 - U[j][im]) );
                
                //-- local increase
                dC[q] = dt * ( Iy + Ix );
                
                //-- chemistry
                reduce(t);
                
                //-- store increase
                I[j][i] = dC[q];
            }
        }
    }

}

void Cell:: update_fields(double t)
{
    for( unit_t j=vtop.y;j>0;--j)
    {
        // core
        for( unit_t i=vtop.x;i>0;--i)
        {
            size_t q = 0;
            for( Specs::iterator p = specs.begin(); p != specs.end(); ++p )
            {
                ++q;
                SpeciesData  &sd = **p;
                Array        &U  = *(sd.U);
                const Array  &I  = *(sd.I);
                C[q]     =  U[j][i] +shrink * I[j][i];
            }
            normalize(t);
            q = 0;
            for( Specs::iterator p = specs.begin(); p != specs.end(); ++p )
            {
                ++q;
                SpeciesData  &sd = **p;
                Array        &U  = *(sd.U);
                U[j][i] = C[q];
            }
        }
        
        // sides @left
        {
            size_t q = 0;
            for( Specs::iterator p = specs.begin(); p != specs.end(); ++p )
            {
                ++q;
                SpeciesData  &sd = **p;
                const Array  &U  = *(sd.U);
                C[q] = (4*U[j][1]-U[j][2])/3.0;
            }
            normalize(t);
            q = 0;
            for( Specs::iterator p = specs.begin(); p != specs.end(); ++p )
            {
                ++q;
                SpeciesData  &sd = **p;
                Array        &U  = *(sd.U);
                U[j][0] = C[q];
            }
        }
        
        // sides @right
        const unit_t im2 = volumes.x - 2;
        const unit_t im1 = volumes.x - 1;
        const unit_t im0 = volumes.x;
        {
            size_t q = 0;
            for( Specs::iterator p = specs.begin(); p != specs.end(); ++p )
            {
                ++q;
                SpeciesData  &sd = **p;
                const Array  &U  = *(sd.U);
                C[q] = (4*U[j][im1]-U[j][im2])/3.0;
            }
            normalize(t);
            q = 0;
            for( Specs::iterator p = specs.begin(); p != specs.end(); ++p )
            {
                ++q;
                SpeciesData  &sd = **p;
                Array        &U  = *(sd.U);
                U[j][im0] = C[q];
            }
        }
        
    }
}

