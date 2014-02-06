#include "yocto/fs/local-fs.hpp"
#include "yocto/math/types.hpp"
#include "yocto/string/tokenizer.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/ptr/auto.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/sort/quick.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/exception.hpp"
#include "yocto/ios/ocstream.hpp"
#include <cmath>

using namespace yocto;
using namespace yocto;

static inline bool is_sep( char C ) throw()
{
    return C==' ' || C == '\t';
}

int main(int argc, char *argv[] )
{
    try
    {
        // Find all Cgel*.log
        vector<double> C;
        vector<double> D;
        vector<string> words;
        
        vfs &fs = local_fs::instance();
        const string dirname = ".";
        auto_ptr<vfs::scanner> scan( fs.new_scanner( dirname ) );
        for( const vfs::entry *ep = scan->next(); ep; ep = scan->next() )
        {
            //std::cerr << '<' << ep->attr << '>' << ep->base_name << std::endl;
            if( ep->has_extension("log") )
            {
                string fn = vfs::get_base_name(ep->base_name);
                fn.trim(8);
                fn.skip(4);
                const double cc = strconv::to<double>(fn,"C");
                std::cerr << "C=" << cc << std::endl;
                C.push_back(cc);
                
                ios::icstream fp(ep->path);
                string line;
                if( !fp.read_line(line) < 0 )
                    throw exception("couldn't open %s", ep->base_name);
                
                words.free();
                tokenizer::split(words, line, is_sep);
                if(words.size() < 3 )
                    throw exception("missing values...");
                
                std::cerr << "D=" << words[3] << std::endl;
                
                D.push_back( strconv::to<double>(words[3],"D") );
            }
        }
        
        co_qsort(C, D);
        
        vector<double> pC;
        vector<double> ratio;
        
        {
            ios::ocstream fp("../Dgel.out",false);
            ios::ocstream fp1("../Dgel1.out",false);

            fp("#C D ratio\n");
            fp1("#pC ratio\n");
            for(size_t i=1;i<=C.size();++i)
            {
                fp("%.4e %.10e %.5e\n", C[i], D[i], D[i]/D[1]);
                if(C[i]>0)
                {
                    const double xx = -log10(C[i]);
                    const double yy = D[i]/D[1];
                    fp1("%.5e %.5e\n", -log10(C[i]), D[i]/D[1]);
                    
                    pC.push_back(xx);
                    ratio.push_back(yy);
                }
            }
        }
        
        

        return 0;
    }
    catch(const exception &e )
    {
        std::cerr << e.what() << std::endl;
        std::cerr << e.when() << std::endl;
    }
    catch(...)
    {
    }
    return -1;
}