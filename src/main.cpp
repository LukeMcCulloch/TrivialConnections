#include <iostream>
using namespace std;

#include "Viewer.h"
#include "Direction.h"
using namespace DDG;

int main( int argc, char** argv )
{
   if( argc != 2 )
   {
      cerr << "usage: " << argv[0] << " in.obj" << endl;
      return 1;
   }

   
   std::cout << "--------- ---------" << std::endl;
   std::cout << "Creating viewer" << std::endl;
   Viewer viewer;
   
   std::cout << "--------- ---------" << std::endl;
   std::cout << "read tri soup?" << std::endl;
   viewer.mesh.read( argv[1] );
   
   std::cout << "--------- ---------" << std::endl;
   std::cout << "init mesh" << std::endl;
   viewer.mesh.init();
   
   std::cout << "--------- ---------" << std::endl;
   std::cout << "create direction field object" << std::endl;
   DirectionField field;
   
   std::cout << "--------- ---------" << std::endl;
   std::cout << "generate direction field" << std::endl;
   field.generate( viewer.mesh, viewer.angle );
   
   std::cout << "--------- ---------" << std::endl;
   std::cout << "viewer init" << std::endl;
   viewer.init();

   return 0;
}

