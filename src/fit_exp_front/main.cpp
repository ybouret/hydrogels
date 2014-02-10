#include "yocto/fs/vfs.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/exception.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/math/fit/lsf.hpp"
#include "yocto/code/utils.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/stat/descr.hpp"
#include "yocto/sort/quick.hpp"
#include "yocto/math/fcn/functions.hpp"
#include "yocto/ptr/shared.hpp"

#include <iostream>

using namespace yocto;
using namespace math;


class Front
{
public:
    explicit  Front( const string &fn ) :
    t(),
    x(),
    x2(),
    z2()
    {
        ios::icstream fp(fn);
        data_set<double> ds;
        ds.use(3, t);
        ds.use(2, x);
        
        ds.load(fp,NULL,1,0);
    }
    
    vector<double> t,x,x2,z2;
    
    typedef shared_ptr<Front> Ptr;
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Front);
};

int main( int argc, char *argv[] )
{
    const char *progname = vfs::get_base_name(argv[0]);
    try
    {
        if(argc<=1)
        {
            throw exception("usage: %s file1.txt ...",progname);
        }
        
        for(int i=1;i<argc;++i)
        {
            Front::Ptr pFront( new Front( argv[i] ) );
            
        }
        
        
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
