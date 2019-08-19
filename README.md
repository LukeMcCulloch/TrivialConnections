# TrivialConnections
Discrete Differential Geometry Exterior Calculus:  Keenan Crane's Trivial Connections Paper


### converting the code from mac to linux:

 * suitesparse Changed: <br/>
//#include <SuiteSparseQR.hpp> <br/>
//#include <umfpack.h>

#include <cholmod.h>

to be <br/>
#include <suitesparse/SuiteSparseQR.hpp> <br/>
#include <suitesparse/umfpack.h>


#include <suitesparse/cholmod.h>

 * Also (probably because I dropped -lmetis?)
 I needed to add "-lumfpack -lamd"
 to the linker for umf libs.

 * GL changed

#include <GLUT/glut.h>

to be

#include <GL/glut.h>

 * Added glew include and "-lGLEW " link

#include <GL/glew.h>

to Viewer.h and Shader.h