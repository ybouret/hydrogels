#include "workspace.hpp"
#include "yocto/spade/variables.hpp"

Workspace:: ~Workspace() throw() {}

Workspace:: Workspace( const Parameters &param, library &l ) :
WorkspaceBase( param.vlayout, param.fields , param.noghost ),
X( mesh.X() ),
handles(l.size())
{
    
    //--------------------------------------------------------------------------
    // construct auxiliary fields
    //--------------------------------------------------------------------------
    const unit_t den = X.upper - X.lower;
    for(unit_t i=X.lower;i<=X.upper;++i)
    {
        mesh.X()[i] = ( (i-X.lower) * param.length ) / den;
    }
    mesh.X()[X.upper] = param.length;
    std::cerr << "X=" << X << std::endl;
    
    //--------------------------------------------------------------------------
    // link arrays to species
    //--------------------------------------------------------------------------
    variables var;
    variables dvar;
    for( library::iterator i=l.begin();i != l.end(); ++i )
    {
        species     &sp = **i;
        SpeciesData &sd = sp.get<SpeciesData>();
        const string spf = sp.name + "F";
        const string spi = sp.name + "I";
        Array       &U  = (*this)[ sp.name ].as<Array>();
        Array       &F  = (*this)[ spf     ].as<Array>();
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

void Workspace:: loadC( array<double> &C, unit_t x ) const
{
    assert( C.size() >= handles.size() );
    load<double>(C, handles, x-X.lower);
}


void Workspace:: load_dC( array<double> &dC, unit_t x ) const
{
    assert( dC.size() >= handles.size() );
    load<double>(dC, handles_dC, x-X.lower);
}


void Workspace:: saveC( const array<double> &C, unit_t x )
{
    save<double>(handles, C, x-X.lower);
}

void Workspace:: save_dC( const array<double> &dC, unit_t x )
{
    save<double>(handles_dC, dC, x-X.lower);
}


