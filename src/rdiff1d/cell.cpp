#include "cell.hpp"

Cell:: ~Cell() throw()
{
}

Cell:: Cell( lua_State *L ) :
Collection(L),
ChemSys(L,*this),
Params(L,*this),
Workspace(vlayout,F,noghost),
ini_side(L,"ini_side", *this),
ini_core(L,"ini_core", *this),
on_side(*this),
in_core(*this),
X( mesh.X() ),
pConc()
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
        D.push_back( sp.data.as<SpeciesData>().D );
    }
    
}

