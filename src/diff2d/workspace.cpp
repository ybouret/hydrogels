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
    const unit_t den = X.upper - X.lower;
    for(unit_t i=X.lower;i<=X.upper;++i)
    {
        mesh.X()[i] = ( (i-X.lower) * param.length.x ) / den;
    }
    mesh.X()[X.upper] = param.length.x;
    std::cerr << "X=" << X << std::endl;
    std::cerr << "Y=" << X << std::endl;

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
