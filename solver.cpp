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
// Function to set the grid-related properties -> Sankar: extend?
template <class mType, class dType>
void solver<mType, dType>::setMesh( ) { // @Sankar, maybe play here with the input
    grid.setupGrid();
    N = grid.noGridPoints;
    h = grid.gridSpacing;
}

// -------------------------------------------------------------------------------------------------
// Function to apply boundary conditions - needs to be extended -> Sankar
template <class mType, class dType>
void solver<mType, dType>::applyBC( Eigen::MatrixBase<mType> &inVec ) {

    input_parser inputParserObj;

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

    // On the bottom
    std::string bottom = inputParserObj.getBCBottom();
    parser.parse(byteCode,bottom,varnames);

    for(auto& index: grid.bdryNodeList.bottom)
    {
        byteCode.var[0] = grid.getCooX(index);
        byteCode.var[1] = grid.getCooY(index);
        inVec[index] = byteCode.run();
    }
    // On the right
    std::string right = inputParserObj.getBCRight();
    parser.parse(byteCode,right,varnames);

    for(auto& index: grid.bdryNodeList.right) {
        byteCode.var[0] = grid.getCooX(index);
        byteCode.var[1] = grid.getCooY(index);
        inVec[index] = byteCode.run();
    }
    // On the top
    std::string top = inputParserObj.getBCTop();
    parser.parse(byteCode,top,varnames);

    for(auto& index: grid.bdryNodeList.top)
    {
        byteCode.var[0] = grid.getCooX(index);
        byteCode.var[1] = grid.getCooY(index);
        inVec[index] = byteCode.run();
    }
    // On the left
    std::string left = inputParserObj.getBCLeft();
    parser.parse(byteCode,left,varnames);

    for(auto& index: grid.bdryNodeList.left)
    {
        byteCode.var[0] = grid.getCooX(index);
        byteCode.var[1] = grid.getCooY(index);
        inVec[index] = byteCode.run();
    }
}

// -------------------------------------------------------------------------------------------------
// Method to build poisson matrix
template <class mType, class dType>
void solver<mType, dType>::buildPoissonMatrix( Eigen::SparseMatrix<dType, Eigen::RowMajor> &poissonMatrix ) { 
    
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
    poissonMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
    //~std::cout << poissonMatrix << std::endl;
}

// -------------------------------------------------------------------------------------------------
// Method to get initial guess for Newton iterations
template <class mType, class dType>
void solver<mType, dType>::getInitGuess( Eigen::MatrixBase<mType> &z ){
    
    // Prepare RHS
    mType b = mType::Zero(N*N);
    applyBC(b);
    
    // Build Poisson-matrix 
    Eigen::SparseMatrix<dType, Eigen::RowMajor> poissonMatrix(N*N, N*N);
    buildPoissonMatrix(poissonMatrix);

    // Solve the Poisson equation using BiCGSTAB 
    Eigen::setNbThreads(numThreads);
    Eigen::BiCGSTAB<Eigen::SparseMatrix<dType, Eigen::RowMajor> > bicgstab;
    bicgstab.setTolerance(TOL_linsolver);
    z = bicgstab.compute(poissonMatrix).solve(b);
    
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

       if (fabs(tmp) > 1e-10)
           outVec[i] = tmp;
       else
           outVec[i] = 0.;
   }

}

// -------------------------------------------------------------------------------------------------
// Action of discrete minSurfOperator on an input vector (STL vector)
template <class mType, class dType> 
template <class dcoType> 
std::vector<dcoType> solver<mType, dType>::minSurfOp_Vector(const std::vector<dcoType> &inVec) {
    std::vector<dcoType> outVec(N*N);
    for(auto& i: grid.innerNodeList) {
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
void solver<mType, dType>::minSurfJac_HardCoded( Eigen::SparseMatrix<dType, Eigen::RowMajor> &Jacobian,
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
         
   // Build sparse matrix from triplets
   Jacobian.setFromTriplets(tripletList.begin(), tripletList.end(),
   [] (const dType &a, const dType &b) { return b; }); 
}

// -------------------------------------------------------------------------------------------------
// Jacobian via AD by hand
// returns application of minSurf-operator + Jacobian
template <class mType, class dType> 
void solver<mType, dType>::minSurfJac_ADByHand( Eigen::SparseMatrix<dType, Eigen::RowMajor> &Jacobian, 
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
         // d(...) = fd(x) // @Chenfei: can't we use the functions defined above??
         dType dx = (inVec[i + 1] - inVec[i - 1]) / (2 * h);
         dType dy = (inVec[i + N] - inVec[i - N]) / (2 * h);
         dType dxx = (inVec[i + 1] - 2 * inVec[i] + inVec[i - 1]) / (h * h);
         dType dyy = (inVec[i + N] - 2 * inVec[i] + inVec[i - N]) / (h * h);
         dType dxy = (inVec[i + 1 + N] + inVec[i - 1 - N] - inVec[i + 1 - N] - inVec[i - 1 + N])
                         / (4 * h * h);
 
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

   // Build sparse matrix from triplets
   Jacobian.setFromTriplets(tripletList.begin(), tripletList.end());

}

// #################################################################################################
// minSurf solving routine

// -------------------------------------------------------------------------------------------------
// Get residual by application of minSurf-operator
template <class mType, class dType>
dType solver<mType, dType>::residual_HardCoded( Eigen::MatrixBase<mType> &resVec,
                                      const Eigen::MatrixBase<mType> &solVec ) {
    // computes residual entries in resVec
    // returns norm of r
    
    // F^h(z^h) = r^h
    // r^h contains the residual, which shall go to zero, in the innerNodeList, and 
    // the boundary information on the bdryNodeList
    minSurfOp(resVec, solVec);
        
    dType res = 0;
    int length = (int) grid.innerNodeList.size();
    #pragma omp parallel for reduction(+:res) if(N>=5*NminParallel)
    for(int index=0; index<length; index++)
        res += pow(resVec[grid.innerNodeList[index]], 2);

    return sqrt(res);
}

// Get residual - AD by hand
template <class mType, class dType>
dType solver<mType, dType>::residual_ADByHand( Eigen::SparseMatrix<dType, Eigen::RowMajor> &Jacobian, 
                                      Eigen::MatrixBase<mType> &resVec,
                                      const Eigen::MatrixBase<mType> &solVec) {
    // computes residual entries in resVec
    // returns norm of r
    
    // F^h(z^h) = r^h
    // r^h contains the residual, which shall go to zero, in the innerNodeList, and 
    // the boundary information on the bdryNodeList
    minSurfJac_ADByHand(Jacobian, resVec, solVec);
    //~std::cout << resVec << std::endl;     
        
    dType res = 0;
    int length = (int) grid.innerNodeList.size();
    #pragma omp parallel for reduction(+:res) if(N>=5*NminParallel)
    for(int index=0; index<length; index++)
        res += pow(resVec[grid.innerNodeList[index]], 2);
    return sqrt(res);
}


// -------------------------------------------------------------------------------------------------
// Functions to run solver
// Hardcoded and AD by hand versions

// Run solver using hardcoded Jacobian
template<class mType, class dType>
void solver<mType, dType>::runSolver_HardCoded( ) {

   mType z = mType::Zero(N*N);
   
 
   std::cout << "Chose to run";
   if (inputParserObj.getPoissonGuess()) {
       std::cout << " with non-trivial initial guess." << std::endl;
       getInitGuess(z); 
   }
   else {
       std::cout << " with trivial initial guess." << std::endl;
       applyBC(z);
   }
  
   dType res;
   mType resVec = mType::Zero(N*N);
   mType dz = mType::Zero(N*N); 
   Eigen::SparseMatrix<dType, Eigen::RowMajor> Jacobian(N*N, N*N);
   
   res = residual_HardCoded(resVec, z);
      
   std::cout << "Starting residual: " << res << std::endl;
   
   dType omega = inputParserObj.getrelaxNewton(); // Relaxation parameter for Newton-Raphson

   unsigned iterationIndex = 0;

   // In case initial guess was not horrifically lucky, run Newton-Raphson
   do { 
       
       // get Jacobian
       minSurfJac_HardCoded(Jacobian, z);
      
       // dz_n = grad[F(z_n)]^-1 * F(z_n)
       // To be played with: preconditioner (MUST), 
       // initial guess (maybe inversion of the Poisson-gradient might also help, but no idea), 
       //     tolerance (MUST).. should not be too high, as our main goal is the result of Newton
       Eigen::BiCGSTAB<Eigen::SparseMatrix<dType, Eigen::RowMajor> > bicgstab;
       bicgstab.setTolerance(TOL_linsolver);
       bicgstab.compute(Jacobian);
       dz = bicgstab.solve(resVec);

       // z_{n+1} = z_n - dz
       z -= omega*dz; 
   
       // get residual and resVec -> F(z_n)
       res = residual_HardCoded(resVec, z);
       
       iterationIndex++;
       if(inputParserObj.getfileFreq() > 0) {
           // Writing the output data
           if( !(iterationIndex%inputParserObj.getfileFreq())) {
               // Writing the vtk files for visualization
               structuredGridWriter<mType,dType>(iterationIndex,z);
               // Writing the residual vs. iteration number into a csv file
               residualWriter<dType>(iterationIndex,res);

           }
       }
       std::cout << "\tIteration: " << iterationIndex  << "\tresidual: " << res << std::endl;

   } while (res > inputParserObj.getTOL_Newton() && 
            iterationIndex < inputParserObj.getmaxIters());

   std::cout << "Stopped after " << iterationIndex << " iterations with a residual of "
             << res << "." << std::endl;
   
   // Writing final output data
   structuredGridWriter<mType,dType>(iterationIndex,z);

   // Writing the residual vs. iteration number of the final solution into a csv file
   residualWriter<dType>(iterationIndex,res);

}

// Run solver using AD by hand
template<class mType, class dType>
void solver<mType, dType>::runSolver_ADByHand( ) {
   mType z = mType::Zero(N*N);
   
   std::cout << "Chose to run";
   if (inputParserObj.getPoissonGuess()) {
       std::cout << " with non-trivial initial guess." << std::endl;
       getInitGuess(z); 
   }
   else {
       std::cout << " with trivial initial guess." << std::endl;
       applyBC(z);
   }
 
   dType res;
   mType resVec = mType::Zero(N*N);
   mType dz = mType::Zero(N*N); 
   Eigen::SparseMatrix<dType, Eigen::RowMajor> Jacobian(N*N, N*N);
   
   // Determine initial residual + Jacobian 
   res = residual_ADByHand(Jacobian, resVec, z);
       
   std::cout << "Starting residual: " << res << std::endl;
   
   dType omega = inputParserObj.getrelaxNewton(); // Relaxation parameter for Newton-Raphson
   
   unsigned iterationIndex = 0;
   // In case initial guess was not horrifically lucky, run Newton-Raphson
   do { 
       // dz_n = grad[F(z_n)]^-1 * F(z_n)
       // To be played with: preconditioner (MUST), 
       // initial guess (maybe inversion of the Poisson-gradient might also help, but no idea), 
       //     tolerance (MUST).. should not be too high, as our main goal is the result of Newton
       Eigen::BiCGSTAB<Eigen::SparseMatrix<dType, Eigen::RowMajor> > bicgstab;
       bicgstab.setTolerance(TOL_linsolver);
       bicgstab.compute(Jacobian);
       dz = bicgstab.solve(resVec);

       // z_{n+1} = z_n - dz
       z -= omega*dz; 
   
       // get residual and resVec -> F(z_n)
       res = residual_ADByHand(Jacobian, resVec, z);
       
       iterationIndex++;
       if(inputParserObj.getfileFreq() > 0) {
           // Writing the output data
           if( !(iterationIndex%inputParserObj.getfileFreq()) ) {
               // Writing the vtk files for visualization
               structuredGridWriter<mType,dType>(iterationIndex,z);
               // Writing the residual vs. iteration number into a csv file
               residualWriter<dType>(iterationIndex,res);
         }
       }
       std::cout << "\tIteration: " << iterationIndex  << "\tresidual: " << res << std::endl;

   } while (res > inputParserObj.getTOL_Newton() && 
            iterationIndex < inputParserObj.getmaxIters());

   std::cout << "Stopped after " << iterationIndex << " iterations with a residual of "
             << res << "." << std::endl;

   // Writing final output data
   structuredGridWriter<mType,dType>(iterationIndex,z);

   // Writing the residual vs. iteration number of the final solution into a csv file
   residualWriter<dType>(iterationIndex,res);

}

// Run solver using AD by dco_c++
template<class mType, class dType>
void solver<mType, dType>::runSolver_ADbyDco( ) {

   mType mz = mType::Zero(N*N);
   std::cout << "Chose to run";
   if (inputParserObj.getPoissonGuess()) {
       std::cout << " with non-trivial initial guess." << std::endl;
       getInitGuess(mz); 
   }
   else {
       std::cout << " with trivial initial guess." << std::endl;
       applyBC(mz);
   }

   std::vector<dType> z(N*N);
   for (int i=0; i<N*N; i++)
       z[i] = mz[i];

   typedef typename dco::gt1s<dType>::type ADtype;
   std::vector<ADtype> yt;
   unsigned iterationIndex = 0;

   // CG solver -- maybe reuse if we use a switch "symmetric"
   /* dType a, b, r2, r2old, res;
   std::vector<dType> dz(N*N, 0.0), dy, p(N*N, 0.0), r;
   do {
       // initialise vector r
       std::vector<ADtype> zt(std::begin(z), std::end(z));
       for (auto& i: grid.innerNodeList)
           dco::derivative(zt)[i] = dz[i]; // init z^(1) = dz
       yt = minSurfOp_Vector(zt); // (y, y^(1)) = F^(1)(x, x^(1));
       for (auto& i: grid.innerNodeList)
           p[i] = -1.0*dco::value(yt)[i]-1.0*dco::derivative(yt)[i];
       r = p;

       // calculate residual
       res = 0;
       for (auto& i: grid.innerNodeList)
           res += pow(dco::value(yt)[i], 2);
       res = std::sqrt(res);
       
       // initiailise |r|^2
       r2 = std::inner_product(r.begin(), r.end(), r.begin(), 0.0);
       r2old = r2;

       // matrix solver (CG)
       while (r2 > 1.e-16) {
           for (auto& i: grid.innerNodeList) dco::derivative(zt)[i] = p[i];
           yt = f(zt);
           std::vector<dType> dy = dco::derivative(yt);
           a = r2old / std::inner_product(p.begin(), p.end(), dy.begin(), 0.0); // a = (r*r)/(p*y^(1))
           for (auto& i: grid.innerNodeList) {
               dz[i] += a*p[i]; // dz += a*p
               r[i] -= a*dy[i]; // r -= a*y^(1)
           }
           r2 = std::inner_product(r.begin(), r.end(), r.begin(), 0.0); // (r*r)
           b = r2/r2old;
           r2old = r2;
           for (auto& i: grid.innerNodeList) p[i] = r[i] + b*p[i]; // p = r+b*p
       }
   
       // update z
       for (auto& i: grid.innerNodeList)
           z[i] += dz[i];

       if( !(iterationIndex%1))
           std::cout << "\tAt iteration " << iterationIndex  << " res is " << res << std::endl;
       iterationIndex++;
   } while (res > 1.e-6 && iterationIndex < 100); */

   // BiCGSTAB solver
   dType rho, rho_new, a, b, w, res, res1;
   std::vector<dType> dz(N*N, 0.0), v(N*N, 0.0), p(N*N, 0.0), r(N*N, 0.0),\
   s(N*N, 0.0), t(N*N, 0.0), r0(N*N, 0.0), z1(N*N, 0.0), y(N*N, 0.0);
   while (iterationIndex < 100) {
       
       std::vector<ADtype> zt(std::begin(z), std::end(z));
       for (auto& i: grid.innerNodeList) dco::derivative(zt)[i] = dz[i];
       yt = minSurfOp_Vector(zt); // (y, y^(1)) = F^(1)(x, x^(1));

       // calculate residual
       res = 0;
       for (auto& i: grid.innerNodeList)
           res += pow(dco::value(yt)[i], 2);
       res = std::sqrt(res);
       if (res < inputParserObj.getTOL_Newton()) break;

       // initialisation
       for (auto& i: grid.innerNodeList)
           p[i] = -1.0*dco::value(yt)[i]-1.0*dco::derivative(yt)[i];
       r = p; r0 = r; std::vector<dType> p(N*N, 0.0); rho = 1; a = 1; w = 1;

       // matrix free solver (BiCGSTAB)
       while (std::inner_product(r.begin(), r.end(), r.begin(), 0.0) > 1.e-16) {
           rho_new = std::inner_product(r0.begin(), r0.end(), r.begin(), 0.0);
           b = (rho_new/rho) * (a/w);
           rho = rho_new;
           for (auto& i: grid.innerNodeList) p[i] = r[i] + b*(p[i]-w*v[i]);
           for (auto& i: grid.innerNodeList) dco::derivative(zt)[i] = p[i];
           yt = minSurfOp_Vector(zt);
           v = dco::derivative(yt);
           a = rho / std::inner_product(r0.begin(), r0.end(), v.begin(), 0.0);
           for (auto& i: grid.innerNodeList) {
               dz[i] += a*p[i];
               s[i] = r[i]-a*v[i];
           }
           // calculate res (intermediate)
           for (auto& i: grid.innerNodeList) z1[i] = z[i] + dz[i];
           y = minSurfOp_Vector(z1);
           res1 = 0;
           for (auto& i: grid.innerNodeList)
               res1 += pow(y[i], 2);
           res1 = std::sqrt(res1);
           if (res1 < inputParserObj.getTOL_Newton()) break;
           
           for (auto& i: grid.innerNodeList) dco::derivative(zt)[i] = s[i];
           yt = minSurfOp_Vector(zt);
           t = dco::derivative(yt);
           w = std::inner_product(t.begin(), t.end(), s.begin(), 0.0) / \
               std::inner_product(t.begin(), t.end(), t.begin(), 0.0);
           for (auto& i: grid.innerNodeList) {
               dz[i] += w*s[i];
               r[i] = s[i]-w*t[i];
           }
       }
   
       // update z
       for (auto& i: grid.innerNodeList)
           z[i] += dz[i];

       iterationIndex++;

       if(inputParserObj.getfileFreq() > 0) {
           // Writing the output data
           if( !(iterationIndex%inputParserObj.getfileFreq())) {
               // Writing the vtk files for visualization
               structuredGridWriter<std::vector<dType>,dType>(iterationIndex, z);
               // Writing the residual vs. iteration number into a csv file
               residualWriter<dType>(iterationIndex,res);

           }
       }
       std::cout << "\tIteration: " << iterationIndex  << "\tresidual: " << res << std::endl;




   }

   std::cout << "Stopped after " << iterationIndex << " iterations with a residual of "
             << res << "." << std::endl;

   // Writing final output data
   structuredGridWriter<std::vector<dType>,dType>(iterationIndex, z);

   // Writing the residual vs. iteration number of the final solution into a csv file
   residualWriter<dType>(iterationIndex,res);


   
}

// Main function to run solver dependent on choice of Jacobian
template<class mType, class dType>
void solver<mType, dType>::runSolver( ) {

    std::cout << "Run minSurf-solver with" ;
    int jacobianOpt = inputParserObj.getjacobianOpt();
    switch(jacobianOpt) {
        case 0:
               std::cout << " hardcoded Jacobian option." << std::endl;
               runSolver_HardCoded();
               break;
        case 1:
               std::cout << " adjoint AD by hand Jacobian option." << std::endl;
               runSolver_ADByHand();
               break;
        case 2:
               std::cout << " DCO tangent Jacobian option with matrix free solver." << std::endl;
               runSolver_ADbyDco();
               break;
    }
} 
