#include "yocto/math/fcn/functions.hpp"
#include "yocto/exception.hpp"
#include "yocto/string/conv.hpp"

using namespace yocto;
using namespace math;


int main(int argc, char *argv[])
{

    try
    {
        for(int i=1;i<argc;++i)
        {
            const double p = strconv::to<double>( argv[i], "p" );
            if(p<=0||p>=2)
                throw exception("invalid iErfc(%g)", p);
            const double x = iqerfc(p);
            std::cerr << "iErfc(" << p << ")=" << x << std::endl;
        }
        
        return 0;
    }
    catch(const exception &e)
    {
        std::cerr << "*** " << e.what() << std::endl;
        std::cerr << "*** " << e.when() << std::endl;
    }
    catch(...)
    {
        std::cerr << "Error" << std::endl;
    }
    return -1;
}