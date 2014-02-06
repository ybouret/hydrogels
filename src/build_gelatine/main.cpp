#include "yocto/fs/local-fs.hpp"
#include "yocto/math/types.hpp"
#include "yocto/string/tokenizer.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/ptr/auto.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/sort/quick.hpp"

using namespace yocto;
using namespace yocto;

int main(int argc, char *argv[] )
{
    try
    {
        // Find all Cgel*.log
        vector<double> C;
        vector<double> D;
        vfs &fs = local_fs::instance();
        const string dirname = ".";
        auto_ptr<vfs::scanner> scan( fs.new_scanner( dirname ) );
        for( const vfs::entry *ep = scan->next(); ep; ep = scan->next() )
        {
            //std::cerr << '<' << ep->attr << '>' << ep->base_name << std::endl;
            if( ep->has_extension("dat") )
            {
                string fn = vfs::get_base_name(ep->base_name);
                fn.trim(4);
                fn.skip(4);
                std::cerr << "---> " << fn << std::endl;
                const double cc = strconv::to<double>(fn,"C");
                std::cerr << "C=" << cc << std::endl;
                C.push_back(cc);
                
                D.push_back(0);
            }
        }
        
        co_qsort(C, D);
        
        std::cerr << "C=" << C << std::endl;
        

        return 0;
    }
    catch(...)
    {
    }
    return -1;
}