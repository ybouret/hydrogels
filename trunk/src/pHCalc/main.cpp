#include "calc.h"
#include "yocto/string/vfs-utils.hpp"
#include "yocto/string/conv.hpp"

using namespace aqueous;

void compute_pH()
{
    
    library lib;
    lib.add("H+",1);
    lib.add("HO-",-1);
    lib.add("Na+",1);
    lib.add("Cl-",-1);
    lib.add("AH",0);
    lib.add("A-",-1);
    
    chemsys cs(lib,1e-7);
    
    equilibrium &water = cs.create( "water", pow(10,-strconv::to_real<double>( input_pKw->value(), "pKw") ) );    
    water.add("H+", 1);
    water.add("HO-",1);
    //std::cerr << water << std::endl;
    
    
    const double Ka = pow(10, -strconv::to_real<double>( input_pKa->value(),"pKa") );
    equilibrium &Ac = cs.create("AH",  Ka); 
    Ac.add("H+",1);
    Ac.add("A-",1);
    Ac.add("AH",-1);
    //std::cerr << Ac << std::endl;
    
    
    cs.build();
    
    initializer ini(lib);
    
    //! add electroneutrality
    ini.electroneutrality();
    
    //! add acid conservation
    {
        constraint &cn = ini.create( strconv::to_real<double>( input_Ca->value(),"Ca") );
        cn.add("AH", 1);
        cn.add("A-", 1);
    }
    
    //! chloride
    {
        constraint &cn = ini.create( strconv::to_real<double>( input_HCl->value(),"HCl") );
        cn.add("Cl-", 1);
    }
    
    //! sodium
    {
        constraint &cn = ini.create( strconv::to_real<double>( input_NaOH->value(),"NaOH") );
        cn.add("Na+", 1);
    }
    
    //! initialize it
    ini(cs,0.0);
    
    const double h  = cs.C[1];
    const double pH = -log10(h);
    
    {
        const string val = vformat("%.4f",pH);
        output_pH->value( val.c_str() );
    }
    
    const double alpha = h/(Ka+h);
    {
        const string val = vformat("%.6f",alpha);
        output_alpha->value( val.c_str() );
    }
    
}

int main( int argc, char *argv[] )
{
	
	const char *progname = _vfs::get_base_name(argv[0]);
	try
	{
		makeCalcWindow()->show();
		return Fl::run();
	}
    catch(const exception &e )
    {
        fl_alert("%s\n%s", e.what(), e.when());
    }
	catch(...)
	{
        fl_alert("unhandled exception in %s",progname);
	}
	return 1;
}
