#ifndef TYPES_INCLUDED
#define TYPES_INCLUDED 1

#include "yocto/spade/array2d.hpp"
#include "yocto/spade/workspace.hpp"
#include "yocto/spade/rmesh.hpp"

using namespace yocto;
using namespace spade;

typedef vertex2D<double>::type           Vertex;
typedef array2D<double>                  Array;
typedef array2D<Vertex>                  VertexArray;
typedef layout2D                         Layout;
typedef rmesh<Layout,double>             Mesh;
typedef fields_setup<Layout>             Fields;
typedef workspace<Layout, rmesh, double> WorkspaceBase;

#endif