#include "yocto/aqueous/lua.hpp"
#include "yocto/lua/lua-state.hpp"
#include "yocto/lua/lua-config.hpp"

#include "yocto/exception.hpp"
#include "yocto/string/vfs-utils.hpp"

#include "yocto/swamp/common.hpp"
#include "yocto/ios/ocstream.hpp"

#include "yocto/code/utils.hpp"

#include "yocto/math/dat/linear.hpp"

#include "yocto/wtime.hpp"

using namespace yocto;
using namespace aqueous;
using namespace swamp;



typedef coord1D         Coord;
typedef array1D<double> Array1D;
typedef layout1D        Layout;
typedef workspace<Layout,double,rmesh> Workspace;


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

class SpeciesData
{
public:
    explicit SpeciesData( double Dvalue ) throw() :
    D( Dvalue ),
    U(0),
    Flux(0),
    Grow(0)
    {
    }
    
    ~SpeciesData() throw() {}
    
    SpeciesData( const SpeciesData &other ) throw() :
    D( other.D ),
    U( other.U ),
    Flux( other.Flux ),
    Grow( other.Grow )
    {
    }
    
    
    const double D; //!< diffusion coefficient
    Array1D *U;     //!< associated field
    Array1D *Flux;  //!< associated flux
    Array1D *Grow;  //!< associated gain
    
private:
    YOCTO_DISABLE_ASSIGN(SpeciesData);
};

class Library : public library
{
public:
    explicit Library( lua_State *L ) : library( sizeof(SpeciesData) )
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
        
        const double D = lua_tonumber(L,-1);
        if( D < 0 ) throw exception("'%s': invalid D value", sp.name.c_str());
        sp.make<SpeciesData,double>(D);
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
        electroneutrality();
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


static inline void workspace2solution( solution &s, const Workspace &w, unit_t i )
{
    component::iterator     q = s.begin();
    for( size_t k = s.size; k>0; --k, ++q )
    {
        component &p = *q;
        p.C = w[ p.name ].as<Array1D>()[ i ];
    }
}

static inline void solution2workspace(  Workspace &w, const solution &s, unit_t i )
{
    component::const_iterator     q = s.begin();
    for( size_t k = s.size; k>0; --k, ++q )
    {
        const component &p = *q;
        w[ p.name ].as<Array1D>()[ i ] = p.C;
    }
}

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
    const double           pH_front;
    
    explicit Parameters( lua_State *L, const library &lib ) :
    __GET_NUMBER(ntop),
    volumes(ntop-1),
    sim_layout(0,ntop),
    sim_ghosts(),
    sim_fields(),
    __GET_NUMBER(gel_length),
    __GET_NUMBER(pH_front)
    {
        if( volumes < 1 )
            throw exception("not enough volumes");
        for( species::iterator i=lib.begin(); i != lib.end(); ++i)
        {
            const species &sp = **i;
            sim_fields.add<Array1D>( sp.name, true);
            sim_fields.add<Array1D>( sp.name+"_flux", false);
            sim_fields.add<Array1D>( sp.name+"_incr", false);
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
    
    vector<SpeciesData> spd;
    Initializer ini_left;    //!< compute what is at left
    Initializer ini_right;   //!< compute what is at right
    Initializer ini_core;    //!< compute what is inside core
    solution    sol_left;    //!< computed
    solution    sol_right;   //!< computed
    solution    sol_core;    //!< computed
    solution    sol_temp;    //!< for chemistry
    array_db   &adb;         //!< this
    Array1D    &x;
    Array1D    &x_half;
    Array1D    &ih;
    Array1D    &ih_half;
    bool        noRightFlux;
    double      weight_coef;
    double      weight_diff;
    double      Ki;
    
    
    explicit Cell( lua_State *L ) :
    Library(L),
    ChemSys(L,*this),
    Parameters(L,*this),
    Workspace( sim_layout, sim_ghosts, sim_fields),
    spd(),
    ini_left(  L, *this, "ini_left"  ),
    ini_right( L, *this, "ini_right" ),
    ini_core(  L, *this, "ini_core"  ), 
    sol_left( *this ),
    sol_right( *this ),
    sol_core( *this ),
    sol_temp( *this ),
    adb( *this ),
    x( mesh.X() ),
    x_half(  adb["x_half" ].as<Array1D>() ),
    ih(      adb["ih"     ].as<Array1D>() ),
    ih_half( adb["ih_half"].as<Array1D>() ),
    noRightFlux(false),
    weight_coef(0),
    weight_diff(0),
    Ki(0)
    {
        
        const ChemSys &cs = *this;
        const equilibrium::ptr *ppEq = cs.search("indic");
        if(!ppEq)
            throw exception("Missing Indicator");
        Ki = (**ppEq).K(0);
        
        std::cerr << "Left Initializer: " << std::endl;
        std::cerr << ini_left << std::endl;
        
        
        std::cerr << "Core Initializer: " << std::endl;
        std::cerr << ini_core << std::endl;
        
        std::cerr << "Right Initializer: " << std::endl;
        std::cerr << ini_right << std::endl;
        
        //----------------------------------------------------------------------
        // initialize chemistry
        //----------------------------------------------------------------------
        ini_left(*this,0.0);
        sol_left.get( C );
        std::cerr << "@left =" << std::endl << sol_left  << std::endl;
        
        
        
        ini_core(*this,0.0);
        sol_core.get( C );
        std::cerr << "@core =" << std::endl << sol_core << std::endl;
        
        ini_right(*this,0.0);
        sol_right.get( C );
        std::cerr << "@right=" << std::endl << sol_right << std::endl;
        
        //----------------------------------------------------------------------
        // special boundary condition
        //----------------------------------------------------------------------
        lua_getglobal(L, "noRightFlux" );
        if( lua_isboolean(L, -1) )
            noRightFlux = lua_toboolean(L, -1) ? true : false ;
        std::cerr << "No Right Flux=" << ( noRightFlux ? "true" : "false") << std::endl;
        
        
        
        //----------------------------------------------------------------------
        // initialize physics
        //----------------------------------------------------------------------
        
        //-- make the grid
        region1D<double>::type r(0,gel_length);
        mesh.regular_map_to(r, *this);
        
        //-- compute auxiliary field x_half
        for( unit_t i=0; i<ntop; ++i )
        {
            x_half[i] = 0.5 * ( x[i] + x[i+1] );
        }
        
        //-- compute volume size
        for( unit_t i=1; i < ntop; ++i )
        {
            ih[i] = 1.0 / (x_half[i]-x_half[i-1]);
        }
        
        // compute distance between mesh vertices
        for(unit_t i=0; i <ntop; ++i )
        {
            ih_half[i]=1.0/(x[i+1] -x[i]);
        }
        
        // initialize
        solution2workspace( *this, sol_left, 0);
        for( unit_t i=1; i < ntop; ++i)
            solution2workspace( *this, sol_core, i);
        solution2workspace( *this, sol_right, ntop);
        
        //----------------------------------------------------------------------
        // hook up species
        //----------------------------------------------------------------------
        Workspace &self = *this;
        for( species::iterator i=lib.begin(); i != lib.end(); ++i)
        {
            const species &sp = **i;
            SpeciesData   &data = sp.get<SpeciesData>();
            data.U    = & self[ sp.name ].as<Array1D>();
            data.Flux = & self[ sp.name + "_flux"].as<Array1D>();
            data.Grow = & self[ sp.name + "_incr"].as<Array1D>();
            spd.push_back(data);
        }
        
        const double dx1 = x[ntop] - x[ntop-1];
        const double dx2 = x[ntop] - x[ntop-2];
        weight_coef = dx2/dx1;
        weight_coef *= weight_coef;
        weight_diff  = weight_coef - 1.0;
        if( noRightFlux )
        {
            sol_right.copy( sol_core );
        }
        
    }
    
    virtual ~Cell() throw()
    {
        
    }
    
    
    //==========================================================================
    // Finite volume for one species
    //==========================================================================
    void computeOneIncrease( double t, double dt, SpeciesData &data )
    {
        
        const double D  =   data.D;
        Array1D     &U  = * data.U;
        Array1D     &F  = * data.Flux;
        Array1D     &G =  * data.Grow;
        
        // flux
        for( unit_t i=0; i < ntop; ++i )
        {
            F[i] = -D * (U[i+1]-U[i]) * ih_half[i];
        }
        
        // - div Flux
        for( unit_t i=1; i < ntop; ++i )
        {
            G[i] = -dt * ih[i] * ( F[i] - F[i-1] );
        }
    }
    
    //==========================================================================
    // reduce increases from composition
    //==========================================================================
    void applyChemistry( double t, unit_t i )
    {
        assert( C.size() == spd.size() );
        assert( dC.size() == spd.size() );
        
        //----------------------------------------------------------------------
        //-- collect raw increases
        //----------------------------------------------------------------------
        for( size_t s=spd.size();s>0;--s )
        {
            SpeciesData &data = spd[s];
            C[s]   = (* data.U)   [i];
            dC[s]  = (* data.Grow)[i];
        }
        
        //----------------------------------------------------------------------
        //-- reduce
        //----------------------------------------------------------------------
        reduce(t);
        
        //----------------------------------------------------------------------
        //-- dispatch raw increases
        //----------------------------------------------------------------------
        for( size_t s=spd.size();s>0;--s )
        {
            SpeciesData &data = spd[s];
            (*data.Grow)[i]   = dC[s];
        }
    }
    
    //==========================================================================
    // compute increases for all species
    //==========================================================================
    void computeIncreases( double t, double dt )
    {
        //std::cerr << "*** Diffusion" << std::endl;
        //----------------------------------------------------------------------
        //-- collect raw increase
        //----------------------------------------------------------------------
        for( size_t s=spd.size();s>0;--s )
        {
            computeOneIncrease(t, dt, spd[s] );
        }
        
        //----------------------------------------------------------------------
        //-- reduce chemistry on each point
        //----------------------------------------------------------------------
        for( size_t i=1; i < ntop; ++i )
        {
            applyChemistry(t, i);
        }
    }
    
    
    //==========================================================================
    // pseudo adaptive step
    //==========================================================================
    void perform( double t, double dt )
    {
        double t_now  = t;  //! where we are
        double dt_rem = dt; //! remaining step
        
        while( dt_rem > 0 )
        {
            computeIncreases( t_now, dt_rem );
            double dt_done = dt_rem;
            double factor  = 1;
            
        CHECK: ;
            for( size_t i=1; i < ntop; ++i )
            {
                for( size_t s=spd.size();s>0;--s )
                {
                    const SpeciesData &data = spd[s];
                    const Array1D &U = *data.U;
                    const Array1D &G = *data.Grow;
                    const double U_i = U[i]; assert( U[i] >= 0 );
                    const double G_i = G[i];
                    if( G_i <= -U_i )
                    {
                        const double tmp = -U_i/(G_i+G_i);
                        if( tmp < factor )
                            factor = tmp;
                    }
                }
            }
            
            if( factor < 1 )
            {
                dt_done *= factor;
                for( size_t i=1; i < ntop; ++i )
                {
                    for( size_t s=spd.size();s>0;--s )
                    {
                        (*spd[s].Grow)[i] *= factor;
                    }
                }
                factor = 1;
                goto CHECK;
            }
            
            
            //------------------------------------------------------------------
            // increase
            //------------------------------------------------------------------
            for( size_t i=1; i < ntop; ++i )
            {
                for( size_t s=spd.size();s>0;--s )
                {
                    SpeciesData &data = spd[s];
                    Array1D       &U = *data.U;
                    const Array1D &G = *data.Grow;
                    U[i] += G[i];
                }
            }
            
            //------------------------------------------------------------------
            // boundary
            //------------------------------------------------------------------
            if( noRightFlux )
            {
                //std::cerr << "weight_coef=" << weight_coef << std::endl;
                //std::cerr << "weight_diff=" << weight_diff << std::endl;
                
                for( size_t s=spd.size();s>0;--s )
                {
                    SpeciesData   &data  = spd[s];
                    Array1D       &U     = *data.U;
                    const double   Utop  = (weight_coef * U[ntop-1] - U[ntop-2])/weight_diff;
                    U[ntop]              = max_of<double>(Utop,0.0);
                    //std::cerr << "Uright #" << s << " = " << U[ntop] << std::endl;
                }
            }
            dt_rem -= dt_done;
            t_now  += dt_done;
            
        }
        
    }
    
    
    bool find_front( double &pos ) const
    {
        const Workspace &W = *this;
        const Array1D   &h = W["H+"].as<Array1D>();
        const c_array<double> arr_x( (double*)&x[0], ntop+1);
        const c_array<double> arr_h( (double*)&h[0], ntop+1);
        vector<double>  front(1,as_capacity);
        
        const double h_front = pow(10.0, -pH_front);
        math::linear_find(h_front, front, arr_x, arr_h, 1);
        if( front.size() )
        {
            pos = front.front();
            return true;
        }
        else
            return false;
    }
    
};

static inline 
void save_profile( const string &filename, const Cell &cell, const string &id )
{
    const Workspace &W = cell;
    const Array1D   &U = W[id].as<Array1D>();
    const Array1D   &x = cell.x;
    ios::ocstream    fp(filename,false);
    for( unit_t i=0; i <= cell.ntop; ++i )
    {
        fp("%g %g\n", x[i], U[i]);
    }
    
}

static inline 
void save_flux( const string &filename, const Cell &cell, const string &id )
{
    const Workspace &W = cell;
    const Array1D   &F = W[id + "_flux" ].as<Array1D>();
    const Array1D   &x = cell.x;
    ios::ocstream    fp(filename,false);
    for( unit_t i=0; i < cell.ntop; ++i )
    {
        fp("%g %g\n", x[i], F[i]);
    }
    
}

static inline 
void save_grow( const string &filename, const Cell &cell, const string &id )
{
    const Workspace &W = cell;
    const Array1D   &G = W[id + "_incr" ].as<Array1D>();
    const Array1D   &x = cell.x;
    ios::ocstream    fp(filename,false);
    for( unit_t i=1; i < cell.ntop; ++i )
    {
        fp("%g %g\n", x[i], G[i]);
    }
    
}

static inline
bool save_front( const string &filename, const Cell &cell, double t )
{
    const bool append = t > 0;
    ios::ocstream fp( filename, append );
    double pos = 0;
    if( cell.find_front(pos) )
    {
        fp("%.15g %.15g\n", t, pos);
        return true;
    }
    else
        return false;
}


#if HAS_FLTK == 1
#include "ui.h"
#include "yocto/auto-ptr.hpp"

static inline void load_pH( const Cell &cell, bool rescale = false )
{
    const Workspace &W = cell;
    FLTK::Curve &crv = Ca->curves["pH"];
    crv.free();
    const Array1D &x = cell.x;
    const Array1D &h = W["H+"].as<Array1D>();
    
    for( unit_t i=0; i <= cell.ntop; ++i )
    {
        crv.push_back( FLTK::Point(x[i],-log10(h[i])) );
    }
    
    Ca->xaxis.set_range(x[0],x[cell.ntop]);
    if( rescale )
        Ca->yaxis.set_range(-log10(h[0]), -log10(h[cell.ntop]) );
    Cmax->value(  -log10(h[cell.ntop]) );
    Ca->redraw();
    
    Scale->xaxis.set_range(x[0],x[cell.ntop]);
    Scale->yaxis.set_range(0, 1);
    
    
    //! fill scale
    FLTK::Points &data = Scale->data;
    data.free();
    for( unit_t i=0; i <= cell.ntop; ++i )
    {
        const double weight = h[i] / (cell.Ki + h[i] );
        //const double weight = h[i] < cell.Ki ? 0 : 1;
        data.push_back( FLTK::Point(x[i],weight ) );
    }
    
    Scale->redraw();
}

static inline void load_alpha( const Cell &cell, bool rescale = false )
{
    const Workspace &W = cell;
    FLTK::Curve &crv = Ca->curves2["alpha"];
    crv.free();
    crv.color = FL_CYAN;
    const Array1D &x = cell.x;
    const Array1D &h = W["H+"].as<Array1D>();
    
    for( unit_t i=0; i <= cell.ntop; ++i )
    {
        const double alpha = h[i] / (cell.Ki+h[i]);
        crv.push_back( FLTK::Point(x[i],alpha) );
    }
    
    Ca->xaxis.set_range(x[0],x[cell.ntop]);
    if( rescale )
        Ca->y2axis.set_range( crv.front().y, crv.back().y );
    Ca->redraw();
    
}



static inline void load_CO2( const Cell &cell )
{
    const Workspace &W = cell;
    FLTK::Curve &crv = Ca->curves2["CO2"];
    crv.free();
    crv.color = FL_RED;
    const Array1D &x = cell.x;
    const Array1D &y = W["CO2"].as<Array1D>();
    
    for( unit_t i=0; i <= cell.ntop; ++i )
    {
        crv.push_back( FLTK::Point(x[i],y[i]) );
    }
    
    
    Ca->xaxis.set_range(x[0],x[cell.ntop]);
    Ca->y2axis.autoscaleY(crv);
    Cmax->value( Ca->y2axis.vmax );
    Ca->redraw();
    
}

#endif

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
        // Create the Cell
        ////////////////////////////////////////////////////////////////////////
        Cell cell(L);
        const bool hasCO2 = cell.lib.search("CO2") != NULL;
        
#if HAS_FLTK
        auto_ptr<Fl_Window> win( makeUI() );
        win->show();
        Tmx->value(0.0);
        Scale->color1 = fl_rgb_color(203, 1,1);
        Scale->color2 = fl_rgb_color(255, 168, 4);
        
        load_pH(cell,true);
        load_alpha(cell,true);
        if( hasCO2 )
            load_CO2(cell);
        Fl::check();
#endif
        
        double Dmax = 0;
        for( library::const_iterator i = cell.lib.begin(); i != cell.lib.end(); ++i )
        {
            const species &sp = **i;
            SpeciesData   &data = sp.get<SpeciesData>();
            if( data.D > Dmax )
                Dmax = data.D;
        }
        std::cerr << "Dmax=" << Dmax << std::endl;
        
        
        double dx_min = fabs(cell.x[cell.ntop] - cell.x[0]);
        for( unit_t i=0; i < cell.ntop; ++i )
        {
            const double tmp = fabs(cell.mesh.dX()[i]);
            if( tmp < dx_min ) dx_min = tmp;
        }
        std::cerr << "dx_min=" << dx_min << std::endl;
        
#define ALPHA_DIFF_MAX (0.1)
        const double D_over_dx2 = Dmax / (dx_min * dx_min );
        const double dt_max     = ALPHA_DIFF_MAX / D_over_dx2;
        const double dt_log     = pow(10.0,floor(log10(dt_max)));
        const double dt_rnd     = floor(dt_max/dt_log);
        
        
        const double dt = dt_rnd * dt_log;
        
        std::cerr << "dt=" << dt << std::endl;
        
        const double dt_sav = 1;
        const size_t every  = max_of<size_t>(1,dt_sav / dt);
        std::cerr << "Saving every " << every << " step" << std::endl;
        
        
        save_front( "front.dat", cell, 0.0 );
        wtime chrono;
        chrono.start();
        double tmx_old = chrono.query();
        std::cerr.flush();
        for( size_t  iter=1; ; ++iter )
        {
            double t_old = (iter-1) * dt;
            double t_now =  iter    * dt;
            cell.perform(t_old,dt);
            if( (iter%every) == 0 )
            {
                fprintf( stderr, "t= %.4g          \r",t_now);fflush(stderr);
                if( !save_front( "front.dat", cell, t_now ) )
                    break;
            }
#if HAS_FLTK
            const double tmx_now = chrono.query();
            if( tmx_now - tmx_old >= 0.1 )
            {
                Tmx->value(t_now);
                load_pH(cell);
                load_alpha(cell);
                if( hasCO2 )
                    load_CO2(cell);
                Fl::check();
                if( !UI_Window->shown() )
                    break;
                tmx_old = tmx_now;
            }
#endif
        }
        fprintf(stderr,"\n\n");
        return 0;
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
	return 1;
}
