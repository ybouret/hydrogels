#include "workspace.hpp"
#include "yocto/spade/variables.hpp"

Workspace:: ~Workspace() throw() {}

Workspace:: Workspace( const Parameters &param, library &l ) :
WorkspaceBase( param.vlayout, param.fields , param.noghost ),
X( mesh.X() ),
Y( mesh.Y() ),
porosity( (*this)[ "porosity"].as<Array>() ),
handles(l.size())
{
    
    //--------------------------------------------------------------------------
    // construct auxiliary fields
    //--------------------------------------------------------------------------
    {
        const unit_t den = X.upper - X.lower;
        for(unit_t i=X.lower;i<=X.upper;++i)
        {
            mesh.X()[i] = ( (i-X.lower) * param.length.x ) / den;
        }
        mesh.X()[X.upper] = param.length.x;
        std::cerr << "X=" << X << std::endl;
    }
    
    {
        const unit_t den = Y.upper - Y.lower;
        for(unit_t i=Y.lower;i<=Y.upper;++i)
        {
            mesh.Y()[i] = ( (i-Y.lower) * param.length.y ) / den;
        }
        mesh.Y()[Y.upper] = param.length.y;
        std::cerr << "Y=" << Y << std::endl;
    }
    
    //--------------------------------------------------------------------------
    // link arrays to species
    //--------------------------------------------------------------------------
    for( library::iterator i=l.begin();i != l.end(); ++i )
    {
        species     &sp = **i;
        SpeciesData &sd = sp.get<SpeciesData>();
        const string spf = sp.name + "F";
        const string spi = sp.name + "I";
        Array       &U  = (*this)[ sp.name ].as<Array>();
        VertexArray &F  = (*this)[ spf     ].as<VertexArray>();
        Array       &I  = (*this)[ spi     ].as<Array>();
        sd.U = &U;
        sd.F = &F;
        sd.I = &I;
        
        var.append(sp.name);
        dvar.append(spi);
    }
    
    //--------------------------------------------------------------------------
    // prepare handles to concentrations
    //--------------------------------------------------------------------------
    query(handles,var);
    query(handles_dC,dvar);
}


void Workspace:: loadC( array<double> &C, const Coord &u ) const
{
    assert( C.size() >= handles.size() );
    load<double>(C, handles, offset_of(u) );
}

void Workspace:: saveC( const array<double> &C, const Coord &u )
{
    assert( C.size() >= handles.size() );
    save<double>(handles, C, offset_of(u) );
}

