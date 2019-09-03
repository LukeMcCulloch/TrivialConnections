/*
 * Trivial Connections on Discrete Surfaces
 * Keenan Crane, Mathieu Desbun and Peter Schroeder
 * SGP 2010 / Computer Graphics Forum
 *
 * A Simplified Algorithm for Simply-Connected Surfaces 
 * Fernando de Goes and Keenan Crane
 *
 * TODO: add soft and hard directional constraints
 */

#ifndef DDG_APPLICATION_H
#define DDG_APPLICATION_H

#include "Mesh.h"
#include "Real.h"
#include "Quaternion.h"
#include "DenseMatrix.h"
#include "SparseMatrix.h"
#include "DiscreteExteriorCalculus.h"
#include "Utility.h"

namespace DDG
{
   class Application
   {
   public:
      bool solveForConnection(Mesh& mesh)
      {
         std::cout << "      ------------------" << std::endl;
         std::cout << "      computing connection" << std::endl;


         bool ok = checkGaussBonnet(mesh);
         if( not ok )
         {
            std::cout << "Gauss-Bonnet thm does not hold" << std::endl;
            return false;
         }
         

         std::cout << "      ------------------" << std::endl;
         std::cout << "      computing trivial holonomy" << std::endl;
         int t0 = clock();
         solveForTrivialHolonomy(mesh);
         int t1 = clock();
         cout << "[trivial] time: " << seconds( t0, t1 ) << "s" << "\n";


         std::cout << "      ------------------" << std::endl;
         std::cout << "      computing nontrivial holonomy" << std::endl;
         t0 = clock();
         solveForNonTrivialHolonomy(mesh);
         t1 = clock();
         cout << "[nontrivial] time: " << seconds( t0, t1 ) << "s" << "\n";
         
         return true;
      }

      double solveForGeodesic(double dt, Mesh& mesh)
      {
         std::cout << "      ------------------" << std::endl;
         std::cout << "      initial condition" << std::endl;
         // initial condiiton
         DenseMatrix<Real> u0;
         int nb = builImpulseSignal(mesh, u0);
         if( nb == 0 ) return 1.0;
         

         std::cout << "      ------------------" << std::endl;
         std::cout << "      Build Hodge * 0Form" << std::endl;
         // DEC
         SparseMatrix<Real> star0;
         HodgeStar0Form<Real>::build( mesh, star0 );
         
         std::cout << "      ------------------" << std::endl;
         std::cout << "      Build Hodge * 1Form" << std::endl;
         SparseMatrix<Real> star1;
         HodgeStar1Form<Real>::build( mesh, star1 );

         std::cout << "      ------------------" << std::endl;
         std::cout << "      Build Exterior Derivative 0Form" << std::endl;
         SparseMatrix<Real> d0;
         ExteriorDerivative0Form<Real>::build( mesh, d0 );
         
         // zero Neumann boundary condition
         SparseMatrix<Real> L = d0.transpose() * star1 * d0;
         
         // make L positive-definite
         L += Real(1.0e-8)*star0;
         
         // heat flow for short interval
         dt *= sqr(mesh.meanEdgeLength());
         SparseMatrix<Real> A = star0 + Real(dt) * L;

         DenseMatrix<Real> u;
         solvePositiveDefinite(A, u, u0);

         // extract geodesic
         computeVectorField(u, mesh);
         
         DenseMatrix<Real> div;
         computeDivergence(mesh, div);

         DenseMatrix<Real> phi;
         solvePositiveDefinite(L, phi, div);

         setMinToZero(phi);
         assignDistance(phi, mesh);         
         return phi.norm();
      }
      
   protected:

   //
   //***********************************************
   // trivial connections
      bool checkGaussBonnet(const Mesh& mesh) const
      {
         // vertex singularity
         int k = 0;
         for(VertexCIter v = mesh.vertices.begin();
             v != mesh.vertices.end();
             v++)
         {
            k += v->singularity;
         }
         
         // generator singularity
         // TODO: include singularity for all generators
         k += mesh.firstGeneratorIndex;
         
         return ( mesh.getEulerCharacteristicNumber() == k );
      }
      
      void solveForTrivialHolonomy(Mesh& mesh)
      {
         // Neumann boundary condition => prescribing geodesic curvature
         // For now, keeping original geodesic curvature
         DenseMatrix<Real> b( mesh.vertices.size() );
         for(VertexIter v = mesh.vertices.begin();
             v != mesh.vertices.end();
             v++)
         {
            double value = 0.0;
            if( not v->onBoundary() )
            {
               value -= ( 2. * M_PI - v->theta() );
               value += 2. * M_PI * v->singularity;
            }
            b( v->index ) = value;
         }
         
         DenseMatrix<Real> u( mesh.vertices.size() );
         if( b.norm() > 1.0e-8 ) backsolvePositiveDefinite( mesh.L, u, b );

         for(VertexIter v = mesh.vertices.begin();
             v != mesh.vertices.end();
             v++)
         {
            v->potential = u( v->index );
         }
      }

      void solveForNonTrivialHolonomy(Mesh& mesh)
      {
         unsigned nb = mesh.numberHarmonicBases();
         if( nb == 0 ) return;
         mesh.harmonicCoefs = std::vector<double>(nb, 0.0);

         DenseMatrix<Real>  b(nb);
         SparseMatrix<Real> H(nb,nb);
         
         int row = 0;
         bool skipBoundaryLoop = true;
         for(unsigned i = 0; i < mesh.generators.size(); ++i)
         {
            const Mesh::Generator& cycle = mesh.generators[i];
            if( skipBoundaryLoop and mesh.isBoundaryGenerator(cycle) )
            {
               skipBoundaryLoop = false;
               continue;
            }
            
            for(unsigned j = 0; j < cycle.size(); ++j)
            {
               HalfEdgeIter he = cycle[j];
               for(unsigned col = 0; col < nb; ++col)
               {
                  H(row,col) += he->harmonicBases[col];
               }
            }

            double value = - mesh.generatorHolonomy( cycle );
            if( row == 0 )
            {
               value += 2.0 * M_PI * mesh.firstGeneratorIndex ;
            }
            b(row) = value;
            row++;
         }
         
         DenseMatrix<Real> x(nb);
         if( b.norm() > 1.0e-8 ) solve(H, x, b);

         for(unsigned i = 0; i < nb; ++i)
            mesh.harmonicCoefs[i] = x(i);
      }

   //
   //***********************************************
   // geodesics in heat
      int builImpulseSignal(const Mesh& mesh, DenseMatrix<Real>& x) const
      {
         int nb = 0;
         x = DenseMatrix<Real>(mesh.vertices.size());
         for( VertexCIter v = mesh.vertices.begin();
             v != mesh.vertices.end();
             v ++ )
         {
            x(v->index) = 0.0;
            if( v->tag )
            {
               x(v->index) = 1.0;
               nb++;
            }
         }
         return nb;
      }
      
      void computeVectorField(const DenseMatrix<Real>& u, Mesh& mesh)
      {
         for( FaceIter f = mesh.faces.begin();
             f != mesh.faces.end();
             f++ )
         {
            if( f->isBoundary() ) continue;
            
            HalfEdgeIter hij = f->he;
            HalfEdgeIter hjk = hij->next;
            HalfEdgeIter hki = hjk->next;

            VertexIter vi = hij->vertex;
            VertexIter vj = hjk->vertex;
            VertexIter vk = hki->vertex;

            double ui = u(vi->index);
            double uj = u(vj->index);
            double uk = u(vk->index);

            Vector eij90 = hij->rotatedEdge();
            Vector ejk90 = hjk->rotatedEdge();
            Vector eki90 = hki->rotatedEdge();

            Vector X = 0.5 * ( ui*ejk90 + uj*eki90 + uk*eij90 ) / f->area();
            f->vector = - X.unit();
         }
      }
      
      void computeDivergence(const Mesh& mesh, DenseMatrix<Real>& div) const
      {
         div = DenseMatrix<Real>(mesh.vertices.size());
         for( VertexCIter v = mesh.vertices.begin();
             v != mesh.vertices.end();
             v ++)
         {
            double sum = 0.0;
            HalfEdgeIter he = v->he;
            do
            {
               if( not he->onBoundary )
               {
                  Vector n = he->next->rotatedEdge();
                  Vector v = he->face->vector;
                  sum += dot( n, v );
               }
               he = he->flip->next;
            }
            while( he != v->he );
            div(v->index) = sum;
         }
      }
      
      void assignDistance(const DenseMatrix<Real>& phi, Mesh& mesh)
      {
         for( VertexIter v = mesh.vertices.begin();
             v != mesh.vertices.end();
             v ++)
         {
            v->distance = phi(v->index);
         }
      }
      
      void setMinToZero(DenseMatrix<Real>& phi) const
      {
         double minValue = 1.0e100;
         for( int i = 0; i < phi.nRows(); ++i )
            for( int j = 0; j < phi.nColumns(); ++j )
               minValue = std::min( minValue, (double) phi(i,j) );

         for( int i = 0; i < phi.nRows(); ++i )
            for( int j = 0; j < phi.nColumns(); ++j )
               phi(i,j) -= minValue;
      }
   };
}

#endif
