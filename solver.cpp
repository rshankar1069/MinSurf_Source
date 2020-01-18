/*  Simulation Sciences Laboratory
 *  WS 2019/20
 *  Chenfei Fan, Praveen Mishra, Sankarasubramanian Ragunathan, Philipp Schleich
 *  Project 1 - "Minimal Surfaces"
 *
 *  Solver class
 */

#include "solver.h"
#include "inputParser.h"
#include "atmsp.h"
#include <fstream>

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

    ATMSP<float> parser;
    ATMSB<float> byteCode;

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

    for(auto& index: grid.bdryNodeList.right)
    {
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
void solver<mType, dType>::buildPoissonMatrix( Eigen::SparseMatrix<dType> &poissonMatrix ) { 
    
    typedef Eigen::Triplet<dType> triplet;
    std::vector<triplet> tripletList;
    tripletList.reserve(5*N);

    // Assemble the FD matrix (should work, validated with matlab)
    // Set inner nodes - 5-pt stencil of FD (there should not be an issue with reaching 
    // the limits of the matrix...)
    for(auto& i: grid.innerNodeList) {
        tripletList.push_back(triplet(i ,i   , 4./(h*h) ));
        tripletList.push_back(triplet(i ,i+1 ,-1./(h*h) ));
        tripletList.push_back(triplet(i ,i-1 ,-1./(h*h) ));
        tripletList.push_back(triplet(i ,i+N ,-1./(h*h) ));
        tripletList.push_back(triplet(i ,i-N ,-1./(h*h) ));
    }

// Set a 1 where Dirichlet BC applies
    // Bottom
    for(auto& i: grid.bdryNodeList.bottom)
        tripletList.push_back(triplet(i, i, 1.));
    // Right
    for(auto& i: grid.bdryNodeList.right)
        tripletList.push_back(triplet(i, i, 1.));
    // Top
    for(auto& i: grid.bdryNodeList.top)
        tripletList.push_back(triplet(i, i, 1.));
    // Left
    for(auto& i: grid.bdryNodeList.left)
        tripletList.push_back(triplet(i, i, 1.));

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
    Eigen::SparseMatrix<dType> poissonMatrix(N*N, N*N);
    buildPoissonMatrix(poissonMatrix);
        
    // Solve the Poisson equation using CG (need to compare a few in the end / fine-tune) 
    Eigen::ConjugateGradient<Eigen::SparseMatrix<dType> > cg;
    z = cg.compute(poissonMatrix).solve(b);
    
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
void solver<mType, dType>::minSurfOperator( Eigen::MatrixBase<mType> &outVec, 
                                                            const Eigen::MatrixBase<mType> &inVec) {
   dType tmp = 0;
   
   for(auto& i: grid.innerNodeList) {
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

       if (fabs(tmp) > .0000000001)
           outVec[i] = tmp;
       else
           outVec[i] = 0.;
   }

}

// -------------------------------------------------------------------------------------------------
// Jacobian by hand
template <class mType, class dType> 
void solver<mType, dType>::minSurfJac_HardCoded( Eigen::SparseMatrix<dType> &Jacobian,
                                                            const Eigen::MatrixBase<mType> &inVec) {
    typedef Eigen::Triplet<dType> triplet;
    std::vector<triplet> tripletList;
    tripletList.reserve(9*N*N);
    
    // Fill Jacobian
    dType dx, dy, dxx, dyy, dxy;
    for(auto& i: grid.innerNodeList) { // is innerNodeList really enough?
        dx = getDx(inVec, i);
        dy = getDy(inVec, i);
        dxx = getDxx(inVec, i);
        dyy = getDyy(inVec, i);
        dxy = getDxy(inVec, i);
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
    }
    
    
    // Build sparse matrix from triplets
    Jacobian.setFromTriplets(tripletList.begin(), tripletList.end(),
    [] (const dType &a, const dType &b) { return b; }); // @Chenfei, I do not really get what this does ^^
}

// -------------------------------------------------------------------------------------------------
// Jacobian by hand
template <class mType, class dType> 
void solver<mType, dType>::minSurfJac_ADByHand( Eigen::SparseMatrix<dType> &Jacobian, 
                                                Eigen::MatrixBase<mType> &outVec, 
                                                const Eigen::MatrixBase<mType> &inVec) {
   const dType h = 1. / ((dType)N);

    typedef Eigen::Triplet<dType> triplet;
    std::vector<triplet> tripletList;
    tripletList.reserve(9 * N * N);

    dType a1_y;
    std::vector<dType> a1_x(9);

    for (auto &i : grid.innerNodeList)
    {
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

        outVec[i] = v1 + v2 + v3;

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
        tripletList.push_back(triplet(i, i - 1 - N, a1_x[0]));
        tripletList.push_back(triplet(i, i     - N, a1_x[1]));
        tripletList.push_back(triplet(i, i + 1 - N, a1_x[2]));
        tripletList.push_back(triplet(i, i - 1,     a1_x[3]));
        tripletList.push_back(triplet(i, i,         a1_x[4]));
        tripletList.push_back(triplet(i, i + 1,     a1_x[5]));
        tripletList.push_back(triplet(i, i - 1 + N, a1_x[6]));
        tripletList.push_back(triplet(i, i     + N, a1_x[7]));
        tripletList.push_back(triplet(i, i + 1 + N, a1_x[8]));
    }
    // Build sparse matrix from triplets
    Jacobian.setFromTriplets(tripletList.begin(), tripletList.end());

}

// #################################################################################################
// minSurf solving routine

// -------------------------------------------------------------------------------------------------
// Get residual by application of minSurf-operator
// Function for hardcoded and ADByHand version

// Get residual - hardcoded
template <class mType, class dType>
dType solver<mType, dType>::residual_HardCoded( Eigen::MatrixBase<mType> &resVec,
                                      const Eigen::MatrixBase<mType> &solVec) {
    // computes residual entries in resVec
    // returns norm of r
    
    // F^h(z^h) = r^h
    // r^h contains the residual, which shall go to zero, in the innerNodeList, and 
    // the boundary information on the bdryNodeList
    minSurfOperator(resVec, solVec);
    //~std::cout << resVec << std::endl;     
        
    dType res = 0;
    for(auto& i: grid.innerNodeList)
        res += pow(resVec[i], 2);

    return sqrt(res);
}

// Get residual - AD by hand
template <class mType, class dType>
dType solver<mType, dType>::residual_ADByHand( Eigen::SparseMatrix<dType> &Jacobian, 
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
    for(auto& i: grid.innerNodeList)
        res += pow(resVec[i], 2);

    return sqrt(res);
}


// -------------------------------------------------------------------------------------------------
// Functions to run solver
// Hardcoded and AD by hand versions

// Run solver using hardcoded Jacobian
template<class mType, class dType>
void solver<mType, dType>::runSolver_HardCoded( ) {
    mType z = mType::Zero(N*N);
    
    getInitGuess(z);
    //~std::cout << z << std::endl;
    //~applyBC(z); // This should be unnecessary if getInitGuess!
    
    dType res;
    mType resVec = mType::Zero(N*N);
    mType dz = mType::Zero(N*N); 
    Eigen::SparseMatrix<dType> Jacobian(N*N, N*N);
    
    res = residual_HardCoded(resVec, z);
        
    std::cout << "Starting residual: " << res << std::endl;
   
    dType omega = 1.0; // Relaxation parameter for Newton-Raphson
    unsigned iterationIndex = 0;
    // In case initial guess was not horrifically lucky, run Newton-Raphson
    do { 
        // get Jacobian
        //~Jacobian.setZero(); // Should be unnecessary
        minSurfJac_HardCoded(Jacobian, z);
        
        // dz_n = grad[F(z_n)]^-1 * F(z_n)
        // To be played with: preconditioner (MUST), 
        // initial guess (maybe inversion of the Poisson-gradient might also help, but no idea), 
        //     tolerance (MUST).. should not be too high, as our main goal is the result of Newton
        Eigen::BiCGSTAB<Eigen::SparseMatrix<dType> > bicgstab;
        bicgstab.compute(Jacobian);
        dz = bicgstab.solve(resVec);

        // z_{n+1} = z_n - dz
        z -= omega*dz; 
    
        // get residual and resVec -> F(z_n)
        res = residual_HardCoded(resVec, z);
        
        iterationIndex++;
        if( !(iterationIndex%1))
            std::cout << "\tAt iteration " << iterationIndex  << " res is " << res << std::endl;
    } while (res > 1.e-6 && iterationIndex < 100);
    std::cout << "Stopped after " << iterationIndex << " iterations with a residual of "
              << res << "." << std::endl;

    std::ofstream file;
    file.open("./results/output.dat",std::ios::out | std::ios::trunc);
    file << z;
    file.close();

}

// Run solver using AD by hand
template<class mType, class dType>
void solver<mType, dType>::runSolver_ADByHand( ) {
    mType z = mType::Zero(N*N);
    
    getInitGuess(z);
    //~std::cout << z << std::endl;
    //~applyBC(z); // This should be unnecessary if getInitGuess!
    
    dType res;
    mType resVec = mType::Zero(N*N);
    mType dz = mType::Zero(N*N); 
    Eigen::SparseMatrix<dType> Jacobian(N*N, N*N);
    
    res = residual_ADByHand(Jacobian, resVec, z);
        
    std::cout << "Starting residual: " << res << std::endl;
   
    dType omega = 1.; // Relaxation parameter for Newton-Raphson
    unsigned iterationIndex = 0;
    // In case initial guess was not horrifically lucky, run Newton-Raphson
    do { 
        // dz_n = grad[F(z_n)]^-1 * F(z_n)
        // To be played with: preconditioner (MUST), 
        // initial guess (maybe inversion of the Poisson-gradient might also help, but no idea), 
        //     tolerance (MUST).. should not be too high, as our main goal is the result of Newton
        Eigen::BiCGSTAB<Eigen::SparseMatrix<dType> > bicgstab;
        bicgstab.compute(Jacobian);
        dz = bicgstab.solve(resVec);

        // z_{n+1} = z_n - dz
        z -= omega*dz; 
    
        // get residual and resVec -> F(z_n)
        res = residual_ADByHand(Jacobian, resVec, z);
        
        iterationIndex++;
        if( !(iterationIndex%1))
            std::cout << "\tAt iteration " << iterationIndex  << " res is " << res << std::endl;
    } while (res > 1.e-6 && iterationIndex < 100);
    std::cout << "Stopped after " << iterationIndex << " iterations with a residual of "
              << res << "." << std::endl;
    //~std::cout << z << std::endl;
}

// -------------------------------------------------------------------------------------------------
// Main solver loop
template<class mType, class dType>
void solver<mType, dType>::runSolver( ) {
    // Determine jacOption from input-file @Sankar
    // ....
    jacOption = 1; // for now...
    
    // Choose how to run solver
    if (jacOption == 0)
        runSolver_HardCoded();
    else if (jacOption == 1)
        runSolver_ADByHand();

}
