# TrivialConnections
Discrete Differential Geometry Exterior Calculus:  Keenan Crane's Trivial Connections Paper


### converting the code 
from mac to linux:

 * suitesparse Changed:
//#include <SuiteSparseQR.hpp>
//#include <umfpack.h>

#include <cholmod.h>

to be 
#include <suitesparse/SuiteSparseQR.hpp>
#include <suitesparse/umfpack.h>


#include <suitesparse/cholmod.h>


 * GL changed

#include <GLUT/glut.h>

to be

#include <GL/glut.h>
