#include "yocto/aqueous/lua.hpp"
#include "yocto/lua/lua-state.hpp"
#include "yocto/lua/lua-config.hpp"

#include "yocto/exception.hpp"
#include "yocto/string/vfs-utils.hpp"

#include "yocto/swamp/common.hpp"

using namespace yocto;
using namespace aqueous;
using namespace swamp;

////////////////////////////////////////////////////////////////////////////////
//
//
//  Chemistry Definitions
//
//
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Library of Species+Diffusion Coefficient
////////////////////////////////////////////////////////////////////////////////

class Library : public library
{
public:
    explicit Library( lua_State *L ) : library( sizeof(double) )
    {
        _lua::species_ctor getD( this, & Library:: loadDiffusionCoefficient );
        _lua::load(L, *this, "species", &getD);
    }
    
    virtual ~Library() throw()
    {
        
    }
    
    
private:
    void loadDiffusionCoefficient( species &sp, lua_State *L )
    {
        if( !lua_isnumber(L,-1) )
        {
            throw exception("'%s': invalid type for D", sp.name.c_str());
        }
        
        double &D = sp.get<double>();
        D = lua_tonumber(L,-1);
        if( D < 0 ) throw exception("'%s': invalid D value", sp.name.c_str());
    }
    YOCTO_DISABLE_COPY_AND_ASSIGN(Library);
};

////////////////////////////////////////////////////////////////////////////////
// Associated chemical system
////////////////////////////////////////////////////////////////////////////////
class ChemSys : public chemsys
{
public:
    explicit ChemSys( lua_State *L, const library &the_library ) :
    chemsys(the_library,1e-5)
    {
        _lua::load(L, *this, "equilibria" );
        build();
    }
    
    
    virtual ~ChemSys() throw()
    {
        
    }
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(ChemSys);
};

////////////////////////////////////////////////////////////////////////////////
// Initializer
////////////////////////////////////////////////////////////////////////////////
class Initializer : public initializer
{
public:
    explicit Initializer( lua_State *L, const library &the_lib, const string &ininame ) :
    initializer(the_lib)
    {
        _lua::load(L, *this, ininame);
    }
    
    virtual ~Initializer() throw()
    {
        
    }
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Initializer);
};

////////////////////////////////////////////////////////////////////////////////
//
//
//  Physics Definition
//
//
////////////////////////////////////////////////////////////////////////////////

typedef coord1D         Coord;
typedef array1D<double> Array1D;
typedef layout1D        Layout;
typedef workspace<Layout,double,rmesh> Workspace;


#define __GET_NUMBER(NAME) NAME( Lua::Config::Get<lua_Number>( L, #NAME ) )

class Parameters 
{
public:
    const unit_t           ntop;
    const unit_t           volumes;
    Layout                 sim_layout;
    ghosts_setup<Coord>    sim_ghosts;
    fields_setup<layout1D> sim_fields;
    const double           gel_length;
    
    explicit Parameters( lua_State *L, const library &lib ) :
    __GET_NUMBER(ntop),
    volumes(ntop-1),
    sim_layout(0,ntop),
    sim_ghosts(),
    sim_fields(),
    __GET_NUMBER(gel_length)
    {
        if( volumes < 1 )
            throw exception("not enough volumes");
        for( species::iterator i=lib.begin(); i != lib.end(); ++i)
        {
            const species &sp = **i;
            sim_fields.add<Array1D>( sp.name, true);
            sim_fields.add<Array1D>( sp.name+"_flux", false);
        }
        sim_fields.add<Array1D>("x_half",false);
        sim_fields.add<Array1D>("ih",false);
        sim_fields.add<Array1D>("h_half",false);
        sim_fields.add<Array1D>("ih_half",false);
    }
    
    virtual ~Parameters() throw()
    {
        
    }
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Parameters);
};


////////////////////////////////////////////////////////////////////////////////
//
//
//  a Cell
//
//
////////////////////////////////////////////////////////////////////////////////
class Cell : 
public Library, 
public ChemSys, 
public Parameters,
public Workspace
{
public:
    
    
    Initializer ini_left;    //!< compute what is at left
    Initializer ini_right;   //!< compute what is at right
    solution    sol_left;    //!< computed
    solution    sol_right;   //!< computed
    array_db   &adb;         //!< this
    Array1D    &x;
    Array1D    &x_half;
    Array1D    &ih;
    Array1D    &ih_half;

    explicit Cell( lua_State *L ) :
    Library(L),
    ChemSys(L,*this),
    Parameters(L,*this),
    Workspace( sim_layout, sim_ghosts, sim_fields),
    ini_left( L, *this, "ini_left" ),
    ini_right( L, *this, "ini_right" ),
    sol_left( *this ),
    sol_right( *this ),
    adb( *this ),
    x( mesh.X() ),
    x_half(  adb["x_half" ].as<Array1D>() ),
    ih(      adb["ih"     ].as<Array1D>() ),
    ih_half( adb["ih_half"].as<Array1D>() )
    {
        //----------------------------------------------------------------------
        // initialize chemistry
        //----------------------------------------------------------------------
        ini_left(*this,0.0);
        sol_left.get( C );
        ini_right(*this,0.0);
        sol_right.get( C );
        
        std::cerr << "@left="  << std::endl << sol_left  << std::endl;
        std::cerr << "@right=" << std::endl << sol_right << std::endl;
        
        //----------------------------------------------------------------------
        // initialize physics
        //----------------------------------------------------------------------
        
        //-- make the grid
        region1D<double>::type r(0,gel_length);
        mesh.regular_map_to(r, *this);
        
        //-- compute auxiliary fields
        for( unit_t i=0; i<ntop; ++i )
        {
            x_half[i] = 0.5 * ( x[i] + x[i+1] );
        }
        
        //-- 
        for( unit_t i=1; i < ntop; ++i )
        {
            ih[i] = 1.0 / (x_half[i]-x_half[i-1]);
        }
        
        for(unit_t i=0; i <ntop; ++i )
        {
            ih_half[i]=1.0/(x[i+1] -x[i]);
        }
    }
    
    virtual ~Cell() throw()
    {
        
    }
};


int main(int argc, char *argv[])
{
    const char *progname = _vfs::get_base_name(argv[0]);
    try 
    {
        
        if( argc < 2 )
            throw exception("Usage: %s config.lua", progname );
        
        ////////////////////////////////////////////////////////////////////////
        // loading parameters
        ////////////////////////////////////////////////////////////////////////
        Lua::State VM;
        lua_State *L = VM();
        Lua::Config::DoFile(L,argv[1]);
        
        ////////////////////////////////////////////////////////////////////////
        // loading chemistry
        ////////////////////////////////////////////////////////////////////////
        
        
        Cell sim(L);
        
        
    }
    catch( const exception &e )
    {
        std::cerr << "in " << progname << std::endl;
        std::cerr << e.what() << std::endl;
        std::cerr << e.when() << std::endl;
    }
    catch(...)
    {
        std::cerr << "unhandled exception in " << progname << std::endl;
    }
	return 0;
}
