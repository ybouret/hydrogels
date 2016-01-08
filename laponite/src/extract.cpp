#include "yocto/program.hpp"
#include "yocto/math/io/data-set.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"

#include "yocto/sequence/vector.hpp"
#include "yocto/fs/local-fs.hpp"

using namespace yocto;
using namespace math;

YOCTO_PROGRAM_START()
{

    if(argc>1)
    {
        const string  filename = argv[1];
        vector<double> labels;
        vector<double> tmx;
        vector<double> pres;
        vector<double> area;

        {
            ios::icstream fp(filename);
            data_set<double> ds;
            ds.use(1,labels);
            ds.use(2,tmx);
            ds.use(3,pres);
            ds.use(4,area);
            ds.load(fp);
        }

        vfs &fs = local_fs::instance();
        const size_t n    = labels.size();
        for(size_t i=1;i<=n;++i)
        {
            const unsigned L       = unsigned(labels[i]);
            const string   new_ext = vformat("b%02u.dat",L);
            string         root    = filename;
            vfs::change_extension(root, new_ext);
            fs.try_remove_file(root);
        }

        for(size_t i=1;i<=n;++i)
        {
            const unsigned L       = unsigned(labels[i]);
            const string   new_ext = vformat("b%02u.dat",L);
            string         root    = filename;
            vfs::change_extension(root, new_ext);
            ios::acstream  fp(root);
            fp("%.15g %.15g %.15g\n", tmx[i], pres[i], area[i]);
        }

    }

}
YOCTO_PROGRAM_END()
