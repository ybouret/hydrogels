#ifndef TYPES_INCLUDED
#define TYPES_INCLUDED 1

#include "yocto/spade/array1d.hpp"
#include "yocto/spade/workspace.hpp"
#include "yocto/spade/rmesh.hpp"

using namespace yocto;
using namespace spade;

typedef array1D<double>                  Array;
typedef layout1D                         Layout;
typedef rmesh<Layout,double>             Mesh;
typedef fields_setup<Layout>             Fields;
typedef workspace<Layout, rmesh, double> WorkspaceBase;

#endif
