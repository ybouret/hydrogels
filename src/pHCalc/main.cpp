#include "calc.h"
#include "yocto/fs/vfs.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/chemical/solution.hpp"

using namespace chemical;

void compute_pH()
{
    
    collection lib;
    lib.add("H+",1);
    lib.add("HO-",-1);
    lib.add("Na+",1);
    lib.add("Cl-",-1);
    lib.add("AH",0);
    lib.add("A-",-1);
    
    equilibria cs;
    
    equilibrium &water = cs.add( "water", pow(10,-strconv::to_real<double>( input_pKw->value(), "pKw") ) );
    water.add( lib["H+"],  1);
    water.add( lib["HO-"], 1);
    //std::cerr << water << std::endl;
    
    
    const double Ka = pow(10, -strconv::to_real<double>( input_pKa->value(),"pKa") );
    equilibrium &Ac = cs.add("AH",  Ka);
    Ac.add( lib["H+"], 1);
    Ac.add( lib["A-"], 1);
    Ac.add( lib["AH"],-1);
    //std::cerr << Ac << std::endl;
    
    
    //cs.build();
    
    boot::loader ini;
    
    //! add electroneutrality
    ini.electroneutrality(lib);
    
    //! add acid conservation
    ini.conserve( lib["AH"], lib["A-"], strconv::to_real<double>( input_Ca->value(),"Ca"));
    
    
    //! chloride
    ini.define( lib["Cl-"], strconv::to_real<double>( input_HCl->value(),"HCl"));
    
    //! sodium
    ini.define( lib["Na+"],strconv::to_real<double>( input_NaOH->value(),"NaOH"));
   
    
    std::cerr << "ini=" << std::endl << ini << std::endl;
    
    //! initialize it
    ini(cs,lib,0.0);
    
    chemical::solution S(lib);
    S.load(cs.C);
    std::cerr << "S=" << S << std::endl;
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
	
	const char *progname = vfs::get_base_name(argv[0]);
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
