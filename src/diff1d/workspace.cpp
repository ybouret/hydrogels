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
    std::cerr << "X=" << X << std::endl;
    
    //--------------------------------------------------------------------------
    // link arrays to species
    //--------------------------------------------------------------------------
    variables var;
    for( library::iterator i=l.begin();i != l.end(); ++i )
    {
        species     &sp = **i;
        SpeciesData &sd = sp.get<SpeciesData>();
        Array       &U  = (*this)[ sp.name ].as<Array>();
        sd.U = &U;
        var.append(sp.name);
    }
    
    //--------------------------------------------------------------------------
    // prepare handles to concentrations
    //--------------------------------------------------------------------------
    query(handles,var);
    
}

void Workspace:: loadC( array<double> &C, unit_t x ) const
{
    assert( C.size() >= handles.size() );
    load<double>(C, handles, x-X.lower);
}


void Workspace:: saveC( const array<double> &C, unit_t x )
{
    save<double>(handles, C, x-X.lower);
}