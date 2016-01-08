#include "yocto/program.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/sequence/vector.hpp"

using namespace yocto;
using namespace math;

YOCTO_PROGRAM_START()
{

    if(argc>1)
    {
        const string  filename = argv[1];
        vector<double> labels;
        vector<double> tmx;
        vector<double> area;

        {
            ios::icstream fp(filename);
            data_set<double> ds;
            ds.use(1,labels);
            ds.use(2,tmx);
            ds.use(3,area);
            ds.load(fp);
        }

    }

}
YOCTO_PROGRAM_END()
