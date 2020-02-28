/*  Simulation Sciences Laboratory
 *  WS 2019/20
 *  Chenfei Fan, Praveen Mishra, Sankarasubramanian Ragunathan, Philipp Schleich
 *  Project 1 - "Minimal Surfaces"
 *
 *  Solver class
 */

#include"solver.h"
#include"inputParser.h"
#include"atmsp.h"
#include"postProcessor.h"
#include<fstream>

// -------------------------------------------------------------------------------------------------
// Function to set the grid-related properties
template <class mType, class dType>
void solver<mType, dType>::setMesh( ) { 
    grid.setupGrid();
    N = grid.noGridPoints;
    h = grid.gridSpacing;
}




// -------------------------------------------------------------------------------------------------
// Function to apply boundary conditions
template <class mType, class dType>
void solver<mType, dType>::applyBC( Eigen::MatrixBase<mType> &inVec ) {

    ATMSP<dType> parser;
    ATMSB<dType> byteCode;

    std::string varnames;

    std::map<std::string,float> consts = inputParserObj.getConsts();
    std::vector<std::string> variables = inputParserObj.getVars();

    for(auto constit=consts.begin(); constit!=consts.end(); constit++)
    {
        parser.addConstant(constit->first,constit->second);
    }
    
    // Storing the variable names in the parser
    
    for(int i=0; i<variables.size(); i++)
    {
        if(i!=variables.size()-1)
            varnames += variables[i]+",";
        else
            varnames += variables[i];
    }

    #pragma omp parallel if(N>=NminParallel)
    {
      #pragma omp sections nowait
      {
        #pragma omp section
        {
          // On the bottom
          std::string bottom = inputParserObj.getBCBottom();
          parser.parse(byteCode,bottom,varnames);

          for(auto& index: grid.bdryNodeList.bottom) {
              byteCode.var[0] = grid.getCooX(index);
              byteCode.var[1] = grid.getCooY(index);
              inVec[index] = byteCode.run();
          }
        }
        #pragma omp section
        {
          // On the right
          std::string right = inputParserObj.getBCRight();
          parser.parse(byteCode,right,varnames);

          for(auto& index: grid.bdryNodeList.right) {
              byteCode.var[0] = grid.getCooX(index);
              byteCode.var[1] = grid.getCooY(index);
              inVec[index] = byteCode.run();
          }
        }
        #pragma omp section
        {
          // On the top
          std::string top = inputParserObj.getBCTop();
          parser.parse(byteCode,top,varnames);

          for(auto& index: grid.bdryNodeList.top) {
              byteCode.var[0] = grid.getCooX(index);
              byteCode.var[1] = grid.getCooY(index);
              inVec[index] = byteCode.run();
          }
        }
        #pragma omp section
        {
          // On the left
          std::string left = inputParserObj.getBCLeft();
          parser.parse(byteCode,left,varnames);

          for(auto& index: grid.bdryNodeList.left) {
              byteCode.var[0] = grid.getCooX(index);
              byteCode.var[1] = grid.getCooY(index);
              inVec[index] = byteCode.run();
          }
        }
      }
    }
}

// -------------------------------------------------------------------------------------------------
// Method to build laplace matrix
template <class mType, class dType>
void solver<mType, dType>::buildLaplaceMatrix( Eigen::SparseMatrix<dType, Eigen::RowMajor> &laplaceMatrix ) { 
    
    typedef Eigen::Triplet<dType> triplet;
    std::vector<triplet> tripletList;
    tripletList.reserve(5*N*N);
 
    int length = (int) grid.innerNodeList.size();
     
    // Assemble the FD matrix (should work, validated with matlab)
    #pragma omp parallel if(N>=5*NminParallel)
    {
       std::vector<triplet> tripletListLoc;
       // Set inner nodes - 5-pt stencil of FD (there should not be an issue with reaching 
       // the limits of the matrix...)
       tripletList.reserve( (int) 6*N*N/numThreads); // 6* here to have a conservative estimate
       #pragma omp for nowait 
       for(int index=0; index < length; index++) {
           int i = grid.innerNodeList[index];
           tripletListLoc.push_back(triplet(i ,i   , 4./(h*h) ));
           tripletListLoc.push_back(triplet(i ,i+1 ,-1./(h*h) ));
           tripletListLoc.push_back(triplet(i ,i-1 ,-1./(h*h) ));
           tripletListLoc.push_back(triplet(i ,i+N ,-1./(h*h) ));
           tripletListLoc.push_back(triplet(i ,i-N ,-1./(h*h) ));
       }
   
       // Set a 1 where Dirichlet BC applies
       #pragma omp sections nowait
       {
         // Bottom
         #pragma omp section
         {
           for(auto& i: grid.bdryNodeList.bottom)
               tripletListLoc.push_back(triplet(i, i, 1.));
         } 
         // Right
         #pragma omp section
         {
           for(auto& i: grid.bdryNodeList.right)
               tripletListLoc.push_back(triplet(i, i, 1.));
         } 
         // Top
         #pragma omp section
         {
           for(auto& i: grid.bdryNodeList.top)
               tripletListLoc.push_back(triplet(i, i, 1.));
         } 
         // Left
         #pragma omp section
         {
           for(auto& i: grid.bdryNodeList.left)
               tripletListLoc.push_back(triplet(i, i, 1.));
         } 
       }
  
        // Combine local triplet list to global one
        #pragma omp critical
        {
          std::move(tripletListLoc.begin(), tripletListLoc.end(),
                   std::back_inserter(tripletList) );
          tripletListLoc.clear();
        }
  
     }
     // Build sparse A from triplets
     laplaceMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
}

// -------------------------------------------------------------------------------------------------
// Method to get initial guess for Newton iterations
// Prescribe average boundary value on the entire grid 
template <class mType, class dType>
void solver<mType, dType>::getInitGuess_Average( Eigen::MatrixBase<mType> &z ){
    
    applyBC(z);
    dType ave=0;
    #pragma omp parallel if(N>=5*NminParallel)
    {
      dType aveLoc=0;
      #pragma omp sections nowait
      {
        #pragma omp section
        {
        for (auto& i: grid.bdryNodeList.bottom) aveLoc += z[i]; 
        }
        #pragma omp section
        {
        for (auto& i: grid.bdryNodeList.right)  aveLoc += z[i]; 
        }
        #pragma omp section
        {
        for (auto& i: grid.bdryNodeList.top)    aveLoc += z[i]; 
        }
        #pragma omp section
        {
        for (auto& i: grid.bdryNodeList.left)   aveLoc += z[i]; 
        }
      }
      #pragma omp atomic
      ave += aveLoc;
    }
    // Compute average over boundary values and set all inner z to average
    ave /= 4*(N-1);
    for (auto& i: grid.innerNodeList) {
      z[i] = ave;
    }

}

// -------------------------------------------------------------------------------------------------
// Method to get Laplace initial guess for Newton iterations
template <class mType, class dType>
void solver<mType, dType>::getInitGuess_Laplace( Eigen::MatrixBase<mType> &z ){
    
    // Prepare RHS
    mType b = mType::Zero(N*N);
    applyBC(b);

    // Build Laplace-matrix 
    Eigen::SparseMatrix<dType, Eigen::RowMajor> laplaceMatrix(N*N, N*N);
    buildLaplaceMatrix(laplaceMatrix);

    // Solve the Laplace equation using BiCGSTAB 
    // using incomplete LU factorization as a preconditioner
    // tuned s.th. outperforms SparseLU decomposition
    Eigen::setNbThreads(numThreads);
    Eigen::BiCGSTAB<Eigen::SparseMatrix<dType, Eigen::RowMajor>,
                    Eigen::IncompleteLUT<dType> >  bicgstab;
    bicgstab.preconditioner().setDroptol(0.001); // preconditioner
    bicgstab.preconditioner().setFillfactor(20); // tuning
    bicgstab.setTolerance(TOL_linsolver);
    z = bicgstab.compute(laplaceMatrix).solve(b);
    
}



// ################################################################################################
// Methods for finite difference stencils

// -------------------------------------------------------------------------------------------------
// Get d(in)/dx by central stencil
template <class mType, class dType>
inline dType solver<mType, dType>::getDx( const Eigen::MatrixBase<mType> &inVec, const int index ) {
    dType out = (inVec[index+1] - inVec[index-1]) / (2*h);
    
    return out;
}

// -------------------------------------------------------------------------------------------------
// Get d(in)/dy by central stencil
template <class mType, class dType>
inline dType solver<mType, dType>::getDy( const Eigen::MatrixBase<mType> &inVec, const int index ) {
    dType out = (inVec[index+N] - inVec[index-N]) / (2*h);
    
    return out;
}

// -------------------------------------------------------------------------------------------------
// Get d^2(in)/dxx by central stencil
template <class mType, class dType>
inline dType solver<mType, dType>::getDxx( const Eigen::MatrixBase<mType> &inVec, const int index ){
    dType out = (inVec[index+1] -2*inVec[index] + inVec[index-1]) / (h*h);
    
    return out;
}

// -------------------------------------------------------------------------------------------------
// Get d^2(in)/dyy by central stencil
template <class mType, class dType>
inline dType solver<mType, dType>::getDyy( const Eigen::MatrixBase<mType> &inVec, const int index ){
    dType out = (inVec[index+N] -2*inVec[index] + inVec[index-N]) / (h*h);
    
    return out;
}

// -------------------------------------------------------------------------------------------------
// Get mixed d^2(in)/dxy by central stencil
template <class mType, class dType>
inline dType solver<mType, dType>::getDxy( const Eigen::MatrixBase<mType> &inVec, const int index ){
    dType out = (inVec[index+1+N] + inVec[index-1-N] - inVec[index+1-N] - inVec[index-1+N])
                                                                                       / (4*h*h);

    return out;
}

// -------------------------------------------------------------------------------------------------
// Action of discrete minSurfOperator on an input vector
template <class mType, class dType> 
void solver<mType, dType>::minSurfOp( Eigen::MatrixBase<mType> &outVec, 
                                                            const Eigen::MatrixBase<mType> &inVec) {
    dType tmp = 0;
 
    int length = (int) grid.innerNodeList.size();
    #pragma omp parallel for firstprivate(tmp) if(N>=5*NminParallel) 
    for(int index=0; index<length; index++) {
        int i = grid.innerNodeList[index];
        // tmp = (1+z_x^2)*z_yy
        tmp = (1 + pow(getDx(inVec, i), 2))
                     * getDyy(inVec, i);
        // tmp -= 2*z_x*z_y*z_xy
        tmp -= 2*  getDx(inVec, i) // z_x
                *  getDy(inVec, i) // z_y
                * getDxy(inVec, i); // z_xy
        // tmp += (1+z_y^2)*z_xx
        tmp += (1 + pow(getDy(inVec, i), 2))
                      * getDxx(inVec, i);

        outVec[i] = tmp;
   }

}

// -------------------------------------------------------------------------------------------------
// Action of discrete minSurfOperator on an input vector (STL vector)
template <class mType, class dType> 
template <class dcoType> 
std::vector<dcoType> solver<mType, dType>::minSurfOp_Vector(const std::vector<dcoType> &inVec) {
    
    std::vector<dcoType> outVec(N*N);
    int length = (int) grid.innerNodeList.size();
    #pragma omp parallel for if(N>=NminParallel)  // tangents are thread-safe
    for(int index=0; index<length; index++) {
        int i = grid.innerNodeList[index];
        // tmp = (1+z_x^2)*z_yy
        outVec[i] = (1 + pow((inVec[i+1] - inVec[i-1]) / (2*h), 2))
                     * (inVec[i+N] -2*inVec[i] + inVec[i-N]) / (h*h) \
        // tmp -= 2*z_x*z_y*z_xy
        - 2*  (inVec[i+1] - inVec[i-1]) / (2*h) // z_x
                *  (inVec[i+N] - inVec[i-N]) / (2*h) // z_y
                * (inVec[i+1+N] + inVec[i-1-N] - inVec[i+1-N] - inVec[i-1+N]) / (4*h*h) // z_xy \
        // tmp += (1+z_y^2)*z_xx
        + (1 + pow((inVec[i+N] - inVec[i-N]) / (2*h), 2))
                      * (inVec[i+1] -2*inVec[i] + inVec[i-1]) / (h*h);
    }
    return outVec;
}

// -------------------------------------------------------------------------------------------------
// Jacobian by hand
template <class mType, class dType> 
void solver<mType, dType>::minSurfJac_Symbolic( Eigen::SparseMatrix<dType, Eigen::RowMajor> &Jacobian,
                                                            const Eigen::MatrixBase<mType> &inVec) {
    typedef Eigen::Triplet<dType> triplet;
    std::vector<triplet> tripletList;
    tripletList.reserve(9*N*N);
    
     // Fill Jacobian
    #pragma omp parallel num_threads(numThreads) if(N>=NminParallel)
    {
      // Divide work among threads -> get numThreads heavy loops and only
      //                              numThreads times communication necessary
      //                              setFromTriplets does not need triplets ordered
      int me = omp_get_thread_num(); 
      int my_start = ( me ) * grid.innerNodeList.size()/numThreads;
      int my_end   = (me+1) * grid.innerNodeList.size()/numThreads;
      if (me == numThreads-1 || N<NminParallel)
          my_end = (int) grid.innerNodeList.size();
      
      // Build triplets locally
      std::vector<triplet> tripletListLoc;
      tripletListLoc.reserve( (int) 10*N*N/numThreads);
      dType dx, dy, dxx, dyy, dxy;
      for(int index = my_start; index < my_end; index++) {
          int i = grid.innerNodeList[index];
          dx = getDx(inVec, i);
          dy = getDy(inVec, i);
          dxx = getDxx(inVec, i);
          dyy = getDyy(inVec, i);
          dxy = getDxy(inVec, i);
          // Build in reduced system
          tripletListLoc.push_back(triplet( i, i, 
                                             -2./(h*h)*(1+pow(dx, 2)) 
                                             -2./(h*h)*(1+pow(dy, 2))
                                             ));
          tripletListLoc.push_back(triplet( i, i+1, 
                                              2./(2*h)*dx*dyy
                                             +1./(h*h)*(1+pow(dy, 2))
                                             -2.*( 1/(2*h) * dy*dxy)
                                             ));
          tripletListLoc.push_back(triplet( i, i-1, 
                                             -2./(2*h)*dx*dyy
                                             +1./(h*h)*(1+pow(dy, 2))
                                             +2.*( 1/(2*h) * dy*dxy)
                                             ));
          tripletListLoc.push_back(triplet( i, i+N,
                                              1./(h*h)*(1+pow(dx, 2))
                                             +2./(2*h)*dy*dxx
                                             -2.*(1/(2*h) * dx*dxy)
                                             ));
          tripletListLoc.push_back(triplet( i, i-N,
                                              1./(h*h)*(1+pow(dx, 2))
                                             -2./(2*h)*dy*dxx
                                             +2.*(1/(2*h) * dx*dxy)
                                             ));
          tripletListLoc.push_back(triplet( i, i+1+N,
                                             -2./(4*h*h)*dx*dy
                                             ));
          tripletListLoc.push_back(triplet( i, i-1+N,
                                              2/(4*h*h)*dx*dy
                                             ));
          tripletListLoc.push_back(triplet( i, i+1-N,
                                              2./(4*h*h)*dx*dy
                                             ));
          tripletListLoc.push_back(triplet( i, i-1-N,
                                        -2./(4*h*h)*dx*dy
                                        ));
          }
          // Combine local triplet list to global one
          #pragma omp critical
          {
            std::move(tripletListLoc.begin(), tripletListLoc.end(),
                      std::back_inserter(tripletList) );
            tripletListLoc.clear();
          }
    }

    // Set a "1" on the diagonal of the boundary nodes to avoid empty rows
    #pragma omp parallel if(N>=NminParallel) 
    {
      #pragma omp sections nowait 
      {
        #pragma omp section
        {
          for(auto& i: grid.bdryNodeList.bottom)
              tripletList.push_back( triplet(i, i, 1) );
        } 
        #pragma omp section
        {
          for(auto& i: grid.bdryNodeList.right)
              tripletList.push_back( triplet(i, i, 1) );
        } 
        #pragma omp section
        {
          for(auto& i: grid.bdryNodeList.top)
              tripletList.push_back( triplet(i, i, 1) );
        } 
        #pragma omp section
        {
          for(auto& i: grid.bdryNodeList.left)
              tripletList.push_back( triplet(i, i, 1) );
        } 
      }
    }
          
    // Build sparse matrix from triplets
    Jacobian.setFromTriplets(tripletList.begin(), tripletList.end(),
    [] (const dType &a, const dType &b) { return b; }); 
}

// -------------------------------------------------------------------------------------------------
// Jacobian via AD by hand
// returns application of minSurf-operator + Jacobian
template <class mType, class dType> 
void solver<mType, dType>::minSurfJac_HandwrittenAdjoint( Eigen::SparseMatrix<dType, 
                                                          Eigen::RowMajor> &Jacobian, 
                                                          Eigen::MatrixBase<mType> &outVec, 
                                                          const Eigen::MatrixBase<mType> &inVec) {
     
    typedef Eigen::Triplet<dType> triplet;
    std::vector<triplet> tripletList;
    tripletList.reserve(9 * N * N);
 
    #pragma omp parallel num_threads(numThreads) if(N>=NminParallel)
    {
      dType a1_y;
      std::vector<dType> a1_x(9);
      // Divide work among threads -> get numThreads heavy loops and only
      //                              numThreads times communication necessary
      //                              setFromTriplets does not need triplets ordered
      int me = omp_get_thread_num(); 
      int my_start = ( me ) * grid.innerNodeList.size()/numThreads;
      int my_end   = (me+1) * grid.innerNodeList.size()/numThreads;
      if (me == numThreads-1 || N<NminParallel)
          my_end = (int) grid.innerNodeList.size();
        
      // Build triplets locally
      std::vector<triplet> tripletListLoc;
      tripletListLoc.reserve( (int) 10*N*N/numThreads);
      // Over local piece of work 
      for(int index = my_start; index < my_end; index++) {
          int i = grid.innerNodeList[index];
          
          // forward
          // d(...) = fd(x)
          dType dx = getDx(inVec, i);
          dType dy = getDy(inVec, i);
          dType dxx = getDxx(inVec, i);
          dType dyy = getDyy(inVec, i);
          dType dxy = getDxy(inVec, i);
  
          // v1 = (1+z_x^2)*z_yy
          dType v1 = (1 + pow(dx, 2)) * dyy;
          // v2 = -2*z_x*z_y*z_xy
          dType v2 = -2 * dx * dy * dxy;
          // v3 = (1+z_y^2)*z_xx
          dType v3 = (1 + pow(dy, 2)) * dxx;
 
          // Forward result - minSurf(inVec)  
          outVec[i] = v1 + v2 + v3; // Index i is only accessed by me <- no atomic necessary
  
          // reverse
          // seed a1_y
          a1_y = 1.;
      
          // reverse of y[i] = v1 + v2 + v3
          dType a1_v1 = a1_y;
          dType a1_v2 = a1_y;
          dType a1_v3 = a1_y;
  
          //reverse of v3 = (1 + pow(dy, 2)) * dxx
          dType a1_dy = dxx * 2 * dy * a1_v3;
          dType a1_dxx = (1 + pow(dy, 2)) * a1_v3;
  
          //reverse v2 = -2 * dx * dy * dxy
          dType a1_dx = -2 * dy * dxy * a1_v2;
          a1_dy += -2 * dx * dxy * a1_v2;
          dType a1_dxy = -2 * dx * dy * a1_v2;
  
          //reverse of v1 = (1 + pow(dx, 2)) * dyy
          a1_dx += dyy * 2 * dx * a1_v1;
          dType a1_dyy = (1 + pow(dx, 2)) * a1_v1;
  
          // reverse of d(...) = fd(x)
          a1_x[5] = a1_dx / (2 * h);
          a1_x[3] = a1_dx / (-2 * h);
          a1_x[7] = a1_dy / (2 * h);
          a1_x[1] = a1_dy / (-2 * h);
          a1_x[5] += a1_dxx / (h * h);
          a1_x[3] += a1_dxx / (h * h);
          a1_x[4] = (-2) * a1_dxx / (h * h);
          a1_x[7] += a1_dyy / (h * h);
          a1_x[1] += a1_dyy / (h * h);
          a1_x[4] += (-2) * a1_dyy / (h * h);
          a1_x[8] = a1_dxy / (4 * h * h);
          a1_x[0] = a1_dxy / (4 * h * h);
          a1_x[2] = a1_dxy / (-4 * h * h);
          a1_x[6] = a1_dxy / (-4 * h * h);
  
          // store derivatives in triplets
          tripletListLoc.push_back(triplet(i, i - 1 - N, a1_x[0]));
          tripletListLoc.push_back(triplet(i, i     - N, a1_x[1]));
          tripletListLoc.push_back(triplet(i, i + 1 - N, a1_x[2]));
          tripletListLoc.push_back(triplet(i, i - 1,     a1_x[3]));
          tripletListLoc.push_back(triplet(i, i,         a1_x[4]));
          tripletListLoc.push_back(triplet(i, i + 1,     a1_x[5]));
          tripletListLoc.push_back(triplet(i, i - 1 + N, a1_x[6]));
          tripletListLoc.push_back(triplet(i, i     + N, a1_x[7]));
          tripletListLoc.push_back(triplet(i, i + 1 + N, a1_x[8]));
      }
     // Combine local triplet list to global one
     #pragma omp critical
     {
       std::move(tripletListLoc.begin(), tripletListLoc.end(),
                 std::back_inserter(tripletList) );
       tripletListLoc.clear();
     }
  
    }  

    // Set a "1" on the diagonal of the boundary nodes to avoid empty rows
    #pragma omp parallel if(N>=NminParallel)
    {
      #pragma omp sections nowait 
      {
       #pragma omp section
       {
         for(auto& i: grid.bdryNodeList.bottom)
              tripletList.push_back( triplet(i, i, 1) );
        } 
        #pragma omp section
        {
          for(auto& i: grid.bdryNodeList.right)
              tripletList.push_back( triplet(i, i, 1) );
        } 
        #pragma omp section
        {
          for(auto& i: grid.bdryNodeList.top)
              tripletList.push_back( triplet(i, i, 1) );
        } 
        #pragma omp section
        {
          for(auto& i: grid.bdryNodeList.left)
              tripletList.push_back( triplet(i, i, 1) );
        } 
      }
    }
 
   // Build sparse matrix from triplets
   Jacobian.setFromTriplets(tripletList.begin(), tripletList.end());

}

// #################################################################################################
// minSurf solving routine

// -------------------------------------------------------------------------------------------------
// Get residual by application of minSurf-operator
template <class mType, class dType>
dType solver<mType, dType>::residual_Symbolic( Eigen::MatrixBase<mType> &resVec,
                                      const Eigen::MatrixBase<mType> &solVec ) {
    // computes residual entries in resVec
    // returns norm of r
    
    // F^h(z^h) = r^h
    // r^h contains the residual, which shall go to zero, in the innerNodeList, and 
    // the boundary information on the bdryNodeList
    minSurfOp(resVec, solVec);
        
    dType res = std::sqrt(innerProduct<mType>(resVec, resVec));
   
    return res;
}

// Get residual - AD by hand (mimicks overloading)
template <class mType, class dType>
dType solver<mType, dType>::residual_HandwrittenAdjoint( Eigen::SparseMatrix<dType, Eigen::RowMajor> &Jacobian, 
                                      Eigen::MatrixBase<mType> &resVec,
                                      const Eigen::MatrixBase<mType> &solVec) {
    // computes residual entries in resVec
    // returns norm of r
    
    // F^h(z^h) = r^h
    // r^h contains the residual, which shall go to zero, in the innerNodeList, and 
    // the boundary information on the bdryNodeList
    minSurfJac_HandwrittenAdjoint(Jacobian, resVec, solVec);
        
    dType res = std::sqrt(innerProduct<mType>(resVec, resVec));

    return res;
}

// Get residual - matrix free 
template <class mType, class dType>
template <class vecType>
dType solver<mType, dType>::residual_matFree( const vecType &resVec) {
    // returns norm of r
    dType res = std::sqrt(innerProduct(resVec, resVec));

    return res;
}



// -------------------------------------------------------------------------------------------------
// Functions to run solver loop
// Options:
// Symbolic differentiation, handwritten adjoint, 
// Matrix-free tangent version using dco_c++ 
 

// Misc functions 
// °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
// Write solver output
template<class mType, class dType>
template<class vecType>
void solver<mType, dType>::writeLoopOutput( const vecType &z, const dType res, const unsigned iteration ) {

   if(inputParserObj.getfileFreq() > 0) {
       // Writing the output data
       if( !(iteration%inputParserObj.getfileFreq()) ) {
           // Writing the vtk files for visualization
           structuredGridWriter<vecType, dType>(iteration, z);
           // Writing the residual vs. iteration number into a csv file
           residualWriter<dType>(iteration,res);

       }
   } 
   std::cout << "\tIteration: " << iteration  << "\tresidual: " << res << std::endl;

}

// Check solver progress
// Returns True, if difference between current and last residual is very small 
// This happens especially using float or for bad tuning parameters
template<class mType, class dType>
bool solver<mType, dType>::checkProgress( const dType res, const dType lastRes ) {

    bool out=false;
    if ( std::fabs(res-lastRes)/res < 1e-2 ) {
        out = true;
    }
    else 
        out = false; 


    return out;

}

// Parallel version of inner product
template<class mType, class dType>
template<class vecType>
dType solver<mType, dType>::innerProduct( const vecType inVec1, const vecType inVec2 ) {
    
    dType out=0;

    // Check if dimensions of inVecs align
    if (inVec1.size() != inVec2.size()) {
        std::cout << "++++++++++ ERROR!! Vectors for computation of inner product have"
                  << " different dimension. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE); 
    }
    int length = (int) inVec1.size();
    #pragma omp parallel for reduction(+:out) if(N>=NminParallel)
    for(int i=0; i<inVec1.size(); i++) {
        out += inVec1[i]*inVec2[i];
    }

    return out;

}


// °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
// Run solver using one of the algorithms that depend on 
// assembling the Jacobian matrix
template<class mType, class dType>
void solver<mType, dType>::runSolver_WithMatrix( Eigen::MatrixBase<mType> &z, int jacobianOpt, 
                                                 dType &res, unsigned &iteration ) {

    res = 0;
    mType resVec = mType::Zero(N*N);
    mType dz  = mType::Zero(N*N);
    Eigen::SparseMatrix<dType, Eigen::RowMajor> Jacobian(N*N, N*N);
    
    // Determine initial residual 
    if (jacobianOpt == 0)
        res = residual_Symbolic(resVec, z);
    else if (jacobianOpt == 1) // also computes Jacobian
        res = residual_HandwrittenAdjoint(Jacobian, resVec, z);
    dType lastRes = res;    
    std::cout << "Starting residual: " << res << std::endl;
    
    dType omega = inputParserObj.getrelaxNewton(); // Relaxation parameter for Newton-Raphson
    
    iteration = 0;
    // In case initial guess was not horrifically lucky, run Newton-Raphson
    do { 

        // if symbolic differentiation, compute Jacobian
        if (jacobianOpt == 0)
            minSurfJac_Symbolic(Jacobian, z);
 
        // dz_n = grad[F(z_n)]^-1 * F(z_n)
        // Solve by BiCGSTAB
        // Why --> testing on the Laplace matrix yielded, that this is 
        // similarly fast as SparseLU, and faster if one uses IncompleteLUT preconditioning
        //
        // Since in a discussion with Klaus we thought, using the full system would not have 
        // any downsides, we pretty much fixated on not to use the reduced system, and 
        // now it is pretty much to late to change this again.
        // BUT: Using the full system yields empty (zero) rows in the Jacobian
        // This is something, BiCGSTAB with Jacobi-Preconditioner (default) can deal with,
        // but LU decomposition cannot. Thus, to increase performance, it would help a 
        // GREAT deal to change to the reduced system! 
        Eigen::BiCGSTAB<Eigen::SparseMatrix<dType, Eigen::RowMajor> > bicgstab;
        bicgstab.setTolerance(TOL_linsolver);
        bicgstab.compute(Jacobian);

        dz = bicgstab.solve(resVec);

        // z_{n+1} = z_n - dz
        z -= omega*dz;
    
        // get residual and resVec -> F(z_n)
        if (jacobianOpt == 0)
            res = residual_Symbolic(resVec, z);
        if (jacobianOpt == 1)
            res = residual_HandwrittenAdjoint(Jacobian, resVec, z);
        iteration++;
        writeLoopOutput<mType>(z, res, iteration);
        // Check if solver stuck
        if (checkProgress(res, lastRes) && iteration>=5) {
            std::cout << "\tSolver is not making any progress. Aborting..." << std::endl;
            break;
        }
        lastRes = res;
 
    } while (res > inputParserObj.getTOL_Newton() && 
             iteration < inputParserObj.getmaxIters());
 

}

// Run solver loop using AD by dco_c++
template<class mType, class dType>
void solver<mType, dType>::runSolver_DcoMatrixFree( mType &mz, dType &res, unsigned &iteration ) {

    std::vector<dType> z(N*N);
    for (int i=0; i<N*N; i++)
        z[i] = mz[i];
 
    typedef typename dco::gt1s<dType>::type ADtype;
    std::vector<ADtype> yt;
    int lengthInnerNodeList = (int) grid.innerNodeList.size();
    bool stopAfterIntermediateLevel = false;                  // Workaround as break statements
                                                              // are not allowed in omp region

    // BiCGSTAB solver 
    // Have non-spd matrices -> CG does not work reliably,
    // Eigen library recommends BiCGSTAB over GMRES
    dType rho, rho_new, a, b, w, res1, lastRes=0;
    std::vector<dType> dz(N*N, 0.0), p(N*N, 0.0), r(N*N, 0.0),\
    s(N*N, 0.0), t(N*N, 0.0), r0(N*N, 0.0), z1(N*N, 0.0), y(N*N, 0.0);
    while (iteration < inputParserObj.getmaxIters()) {
        
        std::vector<ADtype> zt(std::begin(z), std::end(z));
        for (auto& i: grid.innerNodeList) dco::derivative(zt[i]) = dz[i];
        yt = minSurfOp_Vector(zt); // (y, y^(1)) = F^(1)(x, x^(1));
 
        // calculate residual
        res = residual_matFree(dco::value(yt));
        // abort if solver is stuck
        if (checkProgress(res, lastRes) && iteration>=5) {
            std::cout << "\tSolver is not making any progress. Aborting..." << std::endl;
            break;
        }
        lastRes = res; 
        // abort if converged
        if (res < inputParserObj.getTOL_Newton()) 
            break;
 
        // initialisation
        for (auto& i: grid.innerNodeList) {
            r[i] = -1.0*dco::value(yt)[i]-1.0*dco::derivative(yt)[i];
        }
        dType initRes = innerProduct(r, r);
        r0 = r; std::vector<dType> p(N*N, 0.0); std::vector<dType> v(N*N, 0.0);
        rho = 1; a = 1; w = 1;
 
        // matrix free linear iterative solver (BiCGSTAB)
        while (innerProduct(r, r)/initRes > std::pow(inputParserObj.getTOL_linsolver(), 2) || 
                   stopAfterIntermediateLevel) {
            #pragma omp parallel if(N>=NminParallel)
            {
              #pragma omp single
              {
                rho_new = innerProduct(r0, r);
                b = (rho_new/rho) * (a/w);
                rho = rho_new;
              }
              #pragma omp for
              for (int index=0; index < lengthInnerNodeList; index++) {
                    int i = grid.innerNodeList[index];
                    p[i] = r[i] + b*(p[i]-w*v[i]);
                    dco::derivative(zt[i]) = p[i];
              }
              #pragma omp single
              {
                yt = minSurfOp_Vector(zt);
                v = dco::derivative(yt);
                a = rho / innerProduct(r0, v);
              } 
              #pragma omp parallel for if(N>=NminParallel)
              for (int index=0; index < lengthInnerNodeList; index++) {
                    int i = grid.innerNodeList[index];
                    dz[i] += a*p[i];
                    s[i] = r[i]-a*v[i];
              }
              // calculate res (intermediate)
              #pragma omp parallel for if(N>=NminParallel)
              for (int index=0; index < lengthInnerNodeList; index++) {
                    int i = grid.innerNodeList[index];
                    z1[i] = z[i] + dz[i];
              }
              #pragma omp single
              {
                y = minSurfOp_Vector(z1);
                res1 = residual_matFree(y);
                // If solver aborts on intermediate level, final output will be written in runSolver-function and is not lost
                if (res1 < inputParserObj.getTOL_Newton()) {
                    stopAfterIntermediateLevel = true;
                }
                
                for (auto& i: grid.innerNodeList) dco::derivative(zt[i]) = s[i];
                yt = minSurfOp_Vector(zt);
                t = dco::derivative(yt);
                w = innerProduct(t, s) / \
                    innerProduct(t, t);
              }
              #pragma omp parallel for if(N>=NminParallel)
              for (int index=0; index < lengthInnerNodeList; index++) {
                    int i = grid.innerNodeList[index];
                    dz[i] += w*s[i];
                    r[i] = s[i]-w*t[i];
              }
            }
        }
    
        // update z
        for (auto& i: grid.innerNodeList)
            z[i] += dz[i];
 
        iteration++;
        writeLoopOutput(z, res, iteration);
    }
   
}

// Main function to run solver dependent on choice of Jacobian
template<class mType, class dType>
void solver<mType, dType>::runSolver( ) {

    mType z = mType::Zero(N*N);
    
    std::cout << "Chose to run";
    // Choose initial guess
    if (inputParserObj.getInitGuessChoice() == 0) {
        std::cout << " with trivial initial guess." << std::endl;
        applyBC(z);
    }
    else if (inputParserObj.getInitGuessChoice() == 1) {
        std::cout << " with average initial guess." << std::endl;
        getInitGuess_Average(z); 
    }
    else if (inputParserObj.getInitGuessChoice() == 2) {
        std::cout << " with initial guess by solving Laplace's equation." << std::endl;
        getInitGuess_Laplace(z); 
    }
    // Raise error if initial guess is out of range
    else {
        std::cout << "++++++++++ ERROR!! Invalid keyword for initial guess entered. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    }
 
    std::cout << "Run minSurf-solver with" ;
    int jacobianOpt = inputParserObj.getjacobianOpt();
    switch(jacobianOpt) {
        case 0:
               std::cout << " Jacobian by symbolic differentiation." << std::endl;
               break;
        case 1:
               std::cout << " Jacobian by handwritten adjoint AD." << std::endl;
               break;
        case 2:
               std::cout << " DCO tangent Jacobian option with matrix free solver." << std::endl;
               break;
        default: 
               std::cout << "\n\nInvalid option for Jacobian. Try again...";
               exit(EXIT_FAILURE);
               break;
               
    }

    dType res=0; unsigned iteration=0;
    // If using one of the algorithms that assemble a matrix, enter solver loop
    if (jacobianOpt==0 || jacobianOpt==1)
        runSolver_WithMatrix(z, jacobianOpt, res, iteration);
    // If using the matrix-free algorithm, enter the matrix-free solver
    if (jacobianOpt==2)
        runSolver_DcoMatrixFree(z, res, iteration);

    std::cout << "Stopped after " << iteration << " iterations with a residual of "
              << res << "." << std::endl;
 

    // Writing final output data
    structuredGridWriter<mType,dType>(iteration, z);
 
    // Writing the residual vs. iteration number of the final solution into a csv file
    residualWriter<dType>(iteration, res);

    // Set global solution vector (this is needed only for testing)
    solution = z;

} 
