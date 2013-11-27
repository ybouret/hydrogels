#ifndef PARAMS_INCLUDED
#define PARAMS_INCLUDED 1

#include "collection.hpp"
#include "yocto/spade/workspace.hpp"
#include "yocto/spade/in1d.hpp"
#include "yocto/spade/array1d.hpp"

using namespace spade;

typedef coord1D  Coord;
typedef layout1D Layout; 
typedef field_info<Layout>   FieldInfo;
typedef fields_setup<Layout> FieldsSetup;
typedef array1D<double>      Array1D;

class SpeciesParam : public object, public counted
{
public:
    const string name;
    const string flux;
    const string incr;
 
    explicit SpeciesParam( const species &sp );
    virtual ~SpeciesParam() throw();
    const string &key() const throw();
    
    typedef intr_ptr<string,SpeciesParam> Pointer;
    
    typedef set<string,SpeciesParam::Pointer> DB;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(SpeciesParam);
};


class Params : public SpeciesParam::DB
{
public:
    explicit Params( lua_State *L, const collection &lib );
    virtual ~Params() throw();
    
    const size_t       volumes;
    const Layout       vlayout;
    const ghosts_setup noghost;
    const double       length;
    
protected:
    FieldsSetup F;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Params);
};

#endif
