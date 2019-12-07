#include<iostream>
#include<stdio.h>
#include<vector>
#include<valarray>
#include<cmath>
#include<fstream>
#include"Eigen/Eigen/Sparse"
#include"Eigen/Eigen/Eigenvalues"
#include"Eigen/Eigen/Core"
#include"Eigen/Eigen/SparseLU"
#include"Eigen/Eigen/SparseCholesky"
#include"Eigen/Eigen/IterativeLinearSolvers"

#define EPS 1.e-12


// #################################################################################################
// Boundary nodes by lexicographical ordering
// To be extended for a) different boundaries
//                    b) MPI?
template <typename listType>
void setBdryNodes( listType &bdryNodeList, const int N ) {
    int i=0;
    
    // Set first and last row
    for(i=0; i<N; i++) {
    // First row
        bdryNodeList.push_back(i);
    // Last row
        bdryNodeList.push_back(N*(N-1) + i);
//         std::cout << "me: " << i << std::endl;
    }
    
    // outer edges
    for(i=N; i<N*(N-1); i+=N) {
        bdryNodeList.push_back(i);
        bdryNodeList.push_back(i+N-1);
    }
    
    // Sort (works for vector) -> why? my hope is that data is better aligned in the following then
    // and it is worth the effort
    std::sort(bdryNodeList.begin(), bdryNodeList.end());

    // Debug output
    for(auto& i: bdryNodeList)
        std::cout << i << std::endl;
}

// Inner nodes by lexicographical ordering
template <typename listType>
void setInnerNodes( listType &innerNodeList, const int N ) {
    
    // Set inner nodes (assuming structured cartesian grid)
    for(int i=1; i<N-1; i++) {
        for(int j=1; j<N-1; j++) {
            innerNodeList.push_back(i+j*N);
        }
    }
    // Sort list (works for vector)
    std::sort(innerNodeList.begin(), innerNodeList.end());
    // Debug output
    std::cout << "innernodes:\n";
    for(auto& it: innerNodeList)
        std::cout << it << std::endl;
        
}

// Function to apply boundary conditions - needs to be extended -> Sankar
template <typename mType, typename dType, typename listType>
void applyBC( const int BCtype, Eigen::MatrixBase<mType> &V , const listType &innerNodeList, 
            listType bdryNodeList, const int N ) {
    // Set boundary values
    int i,j; dType x, y;
    for(auto& node: bdryNodeList) {
        i = node%N;
        x = ((dType) i) / (N-1);
        j = node/N;
        y = ((dType) j) / (N-1);
        
        if( i==0 or i==(N-1)) // lower or upper bdry -> x^2
            V[node] = pow(y-.5,2);
        else if ( j==0 or j==(N-1)) // left or right bdry -> y^2
            V[node] = pow(x-.5,2);
    }
    //~// Set rest to zero
    //~for(auto& i: innerNodeList)
        //~V[i] = 0;
}

// ################################################################################################
template <typename mType, typename dType, typename listType>
void buildPoissonMatrix( Eigen::SparseMatrix<dType> &A, const listType &innerNodeList, 
                         const listType &bdryNodeList, const int N ) { 

    typedef Eigen::Triplet<dType> triplet;
    std::vector<triplet> tripletList;
    tripletList.reserve(5*N);

    // Assemble the FD matrix (should work, validated with matlab)
    // Set inner nodes - 5-pt stencil of FD (there should not be an issue with reaching 
    // the limits of the matrix...)
    for(auto& i: innerNodeList) {
        tripletList.push_back(triplet(i ,i   , 4));
        tripletList.push_back(triplet(i ,i+1 ,-1));
        tripletList.push_back(triplet(i ,i-1 ,-1));
        tripletList.push_back(triplet(i ,i+N ,-1));
        tripletList.push_back(triplet(i ,i-N ,-1));
    }

// Set a 1 where Dirichlet BC applies
    for(auto& i: bdryNodeList)
        tripletList.push_back(triplet(i, i, 1));

    // Build sparse A from triplets
    A.setFromTriplets(tripletList.begin(), tripletList.end());
    //~std::cout << A << std::endl;
}

template <typename mType, typename dType, typename listType>
void getInitGuess( Eigen::MatrixBase<mType> &zE, const Eigen::MatrixBase<mType> &bE, 
                   const listType &innerNodeList, const listType &bdryNodeList, const int N ){

    // Preallocate Poisson-matrix 
    Eigen::SparseMatrix<dType> A(N*N, N*N);
    
    buildPoissonMatrix<mType, dType, listType>(A, innerNodeList, bdryNodeList, N);
        
    // Now solve the Poisson equation using sparse Cholesky factorization
    // for later: solveWithGuess()
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<dType> > solver;
    zE = solver.compute(A).solve(bE);
    
}

// ################################################################################################

// Get d(in)/dx by central stencil
template <typename mType, typename dType>
inline dType getDx( const Eigen::MatrixBase<mType> &inVec, const dType h, const int index ) {
    dType out = (inVec[index+1] - inVec[index-1]) / (2*h);
    
    return out;
}

// Get d(in)/dy by central stencil
template <typename mType, typename dType>
inline dType getDy( const Eigen::MatrixBase<mType> &inVec, const dType h, const int N, 
                    const int index ) {
    dType out = (inVec[index+N] - inVec[index-N]) / (2*h);
    
    return out;
}

// Get d^2(in)/dxx by central stencil
template <typename mType, typename dType>
inline dType getDxx( const Eigen::MatrixBase<mType> &inVec, const dType h, const int index ) {
    dType out = (inVec[index+1] -2*inVec[index] + inVec[index-1]) / (h*h);
    
    return out;
}

// Get d^2(in)/dyy by central stencil
template <typename mType, typename dType>
inline dType getDyy( const Eigen::MatrixBase<mType> &inVec, const dType h, const int N, 
                     const int index ) {
    dType out = (inVec[index+N] -2*inVec[index] + inVec[index-N]) / (h*h);
    
    return out;
}

// Get mixed d^2(in)/dxy by central stencil
template <typename mType, typename dType>
inline dType getDxy( const Eigen::MatrixBase<mType> &inVec, const dType h, const int N, 
                     const int index ) {
    dType out = (inVec[index+1+N] + inVec[index-1-N] - inVec[index+1-N] - inVec[index-1+N])
                                                                                       / (4*h*h);

    return out;
}



// Action of discrete minSurfOperator on an input vector
template <typename mType, typename dType, typename listType> 
void solver<mType, dType>::minSurfOperator( Eigen::MatrixBase<mType> &outVec,
                                                            const Eigen::MatrixBase<mType> &inVec ){
   const dType h = grid.gridSpacing;
   dType tmp = 0;
   
   outVec.setZero();

   for(auto& i: grid.innerNodeList) {
       // tmp = (1+z_x^2)*z_yy
       tmp = (1 + pow(getDx(inVec, h, i), 2))
                    * getDyy(inVec, h, N, i);
       // tmp -= 2*z_x*z_y*z_xy
       tmp -= 2*  getDx(inVec, h,    i) // z_x
               *  getDy(inVec, h, N, i) // z_y
               * getDxy(inVec, h, N, i); // z_xy
       // tmp += (1+z_y^2)*z_xx
       tmp += (1 + pow(getDy(inVec, h, N, i), 2))
                     * getDxx(inVec, h, i);

       if (fabs(tmp) > .0000000001)
           outVec[i] = tmp;
       else
           outVec[i] = 0.;
   }
   
}

// Jacobian by hand
template <typename mType, typename dType, typename listType> 
void solver<mType, dType>::minSurfJacByHand( Eigen::SparseMatrix<dType> &Jacobian,
                                                            const Eigen::MatrixBase<mType> &inVec) {
    typedef Eigen::Triplet<dType> triplet;
    std::vector<triplet> tripletList;
    tripletList.reserve(9*N*N);
    
    dType h = grid.gridSpacing; // might make sense to store that "globally" 
    // (like in a grid class, grid.h for instance. Then we could have more stuff, such as grid.nNodes, 
    // grid.nInnerNodes, grid.innerNodeList etc, which would not be needed to always handed over)
    // Fill Jacobian
    dType dx, dy, dxx, dyy, dxy;
    for(auto& i: innerNodeList) { // is innerNodeList really enough?
        dx = getDx(inVec, h,    i);
        dy = getDy(inVec, h, N, i);
        dxx = getDxx(inVec, h,    i);
        dyy = getDyy(inVec, h, N, i);
        dxy = getDxy(inVec, h, N, i);
        tripletList.push_back(triplet( i ,i,
                                      -2./(h*h)*(1+pow(dx, 2)) 
                                      -2./(h*h)*(1+pow(dy, 2))
                                      ));
        tripletList.push_back(triplet( i ,i+1, 
                                       2./(2*h)*dx*dyy
                                      +1./(h*h)*(1+pow(dy, 2))
                                      -2.*( 1/(2*h) * dy*dxy)
                                      ));
        tripletList.push_back(triplet( i ,i-1, 
                                      -2./(2*h)*dx*dyy
                                      +1./(h*h)*(1+pow(dy, 2))
                                      +2.*( 1/(2*h) * dy*dxy)
                                      ));
        tripletList.push_back(triplet( i, i+N,
                                       1./(h*h)*(1+pow(dx, 2))
                                      +2./(2*h)*dy*dxx
                                      -2.*(1/(2*h) * dx*dxy)
                                      ));
        tripletList.push_back(triplet( i, i-N,
                                       1./(h*h)*(1+pow(dx, 2))
                                      -2./(2*h)*dy*dxx
                                      +2.*(1/(2*h) * dx*dxy)
                                      ));
        tripletList.push_back(triplet( i ,i+1+N,
                                      -2./(4*h*h)*dx*dy
                                      ));
        tripletList.push_back(triplet( i ,i-1+N,
                                       2/(4*h*h)*dx*dy
                                      ));
        tripletList.push_back(triplet( i ,i+1-N,
                                       2./(4*h*h)*dx*dy
                                      ));
        tripletList.push_back(triplet( i ,i-1-N,
                                      -2./(4*h*h)*dx*dy
                                      ));
        // maybe not useful ...
        /* for(auto& j: bdryNodeList)
            tripletList.push_back(triplet(i, j, 0)); */
    }
    
    
    // Build sparse matrix from triplets
    Jacobian.setFromTriplets(tripletList.begin(), tripletList.end(),
    [] (const dType &a, const dType &b) { return b; }); // @Chenfei, I do not really get what this does ^^
}


// #################################################################################################
// minSurf solving routine

// Get residual by application of minSurf-operator
template <typename mType, typename dType, typename listType>
dType solver<mType, dType>::residual( Eigen::MatrixBase<mType> &resVec,
                                                           const Eigen::MatrixBase<mType> &solVec) {
    // computes residual entries in resVec
    // returns norm of r
    
    // F^h(z^h) = r^h
    // r^h contains the residual, which shall go to zero, in the innerNodeList, and 
    // the boundary information on the bdryNodeList
    resVec.setZero();
    
    minSurfOperator<mType, dType>(resVec, solVec);
    //~std::cout << resVec << std::endl;     
    
    // maybe one could out-source this applyBC, since it is not touched again... 
    // but then, no setting resVec to zero and careful with setZero in minSurfOperator
    //~applyBC<mType, dType, listType>(0, resVec, innerNodeList, bdryNodeList, N); 
    //~std::cout << resVec << std::endl;     
    
    
    dType res = 0;
    for(auto& i: innerNodeList)
        res += pow(resVec[i], 2);

    return sqrt(res);
}


// Main solver loop
template<typename mType, typename dType, typename listType>
void runSolver( const listType &innerNodeList, const listType &bdryNodeList, const int N) {
    mType z = mType::Zero();
    mType b = mType::Zero();
    applyBC<mType, dType, listType>(0, b, innerNodeList, bdryNodeList, N);
    
    getInitGuess<mType, dType, listType>(z, b, innerNodeList, bdryNodeList, N);
    //~std::cout << z << std::endl;
    applyBC<mType, dType, listType>(0, z, innerNodeList, bdryNodeList, N); 
    //~std::cout << z << std::endl;
    
    dType res;
    mType resVec, dz; 
    Eigen::SparseMatrix<dType> Jacobian(N*N, N*N);
    
    res = residual<mType, dType, listType>(resVec, z);
        
    std::cout << "Starting residual: " << res << std::endl;
   
    dType omega = 0.85; // Relaxation parameter for Newton-Raphson, check what might be optimal here
    unsigned iterationIndex = 0;
    // In case initiall guess was not horrifically lucky, run Newton-Raphson
    do { 
        
        // get Jacobian
        Jacobian.setZero();
        minSurfJacByHand(Jacobian, z);
        //~std::cout << Jacobian << std::endl;
        // Test for Eigenvalues of Jacobian - only test purpose, to know whether CG is a good idea or not
        
        // dz_n = grad[F(z_n)]^-1 * F(z_n)
        // To be played with: preconditioner (MUST), 
        // initial guess (maybe inversion of the Poisson-gradient might also help, but no idea), 
        //     tolerance (MUST).. should not be too high, as our main goal is the result of Newton
        Eigen::BiCGSTAB<Eigen::SparseMatrix<dType> > bicgstab;
        bicgstab.compute(Jacobian);
        dz.setZero();
        dz = bicgstab.solve(resVec);
        //~std::cout << dz << std::endl;
        // z_{n+1} = z_n - dz
        z -= omega*dz;   
        applyBC<mType, dType, listType>(z); // is this necessary?
    
        // get residual and resVec -> F(z_n)
        res = residual<mType, dType>(resVec, z);
        
        iterationIndex++;
        if( !(iterationIndex%100))
            std::cout << "\tAt iteration " << iterationIndex  << " res is " << res << std::endl;
    } while (res > 1.e-6 && iterationIndex < 100);
    std::cout << "Stopped after " << iterationIndex << " iterations with a residual of "
              << res << "." << std::endl;
    std::cout << z << std::endl;
}
