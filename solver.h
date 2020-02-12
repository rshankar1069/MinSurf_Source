/*  Simulation Sciences Laboratory
 *  WS 2019/20
 *  Chenfei Fan, Praveen Mishra, Sankrarasubramanian Ragunathan, Philipp Schleich
 *  Project 1 - "Minimal Surfaces"
 *
 *  Solver class -- header file
 */

#include<iostream>
#include<stdio.h>
#include<vector>
#include<numeric>
#include<cmath>
#include"Eigen/Eigen/Core"
#include"Eigen/Eigen/Sparse"
#include"Eigen/Eigen/IterativeLinearSolvers"
#include"inputParser.h"
#include"cartesianGrid.h"
#pragma once

#define DCO_AUTO_SUPPORT
#define DCO_DISABLE_AVX2_WARNING
#include"dco_cpp/include/dco.hpp" 
//#include<dco.hpp> 



template<class mType, class dType> class solver
{
    typedef std::vector<int> listType;
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
        // Constructor
        solver();
        // Destructor
        ~solver();
        
        cartesianGrid<dType, listType> grid;
        input_parser inputParserObj;
        
        void setMesh( );
        void runSolver( );
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    private:
        // Private variables
        // -----------------
        int N; // number of gridpoints in one dimension
        dType h; // grid spacing
        int numThreads; // number of threads for OpenMP multithreading
        int NminParallel; // minimum number of gridpoints for parallel exec
        dType TOL_linsolver; // tolerance of linear iterative solvers

        // Private methods
        // --------------- 
        // Boundary conditions
        void applyBC( Eigen::MatrixBase<mType> &inVec );
        
        // Get initial guess by solving -u''=0 in D, u=g on dD via central FD
        void buildLaplaceMatrix( Eigen::SparseMatrix<dType, Eigen::RowMajor> &laplaceMatrix );
        void getInitGuess_Average( Eigen::MatrixBase<mType> &zE );
        void getInitGuess_Laplace( Eigen::MatrixBase<mType> &zE );
        
        // Helper functions for FD stencils
        inline dType getDx( const Eigen::MatrixBase<mType> &inVec, const int index );
        inline dType getDy( const Eigen::MatrixBase<mType> &inVec, const int index );
        inline dType getDxx( const Eigen::MatrixBase<mType> &inVec, const int index );
        inline dType getDyy( const Eigen::MatrixBase<mType> &inVec, const int index );
        inline dType getDxy( const Eigen::MatrixBase<mType> &inVec, const int index );
        // minSurf-Operator
        void minSurfOp( Eigen::MatrixBase<mType> &outVec,
                        const Eigen::MatrixBase<mType> &inVec );
        // minsurf-Operator (for STL vector)
        template<class dcoType>
        std::vector<dcoType> minSurfOp_Vector(const std::vector<dcoType> &inVec);
        // Jacobian of minSurf-Operator - hardcoded and AD by hand version
        void minSurfJac_Symbolic( Eigen::SparseMatrix<dType, Eigen::RowMajor> &Jacobian,
                               const Eigen::MatrixBase<mType> &inVec);
        void minSurfJac_HandwrittenAdjoint( Eigen::SparseMatrix<dType, Eigen::RowMajor> &Jacobian,
                                 Eigen::MatrixBase<mType> &outVec,
                                 const Eigen::MatrixBase<mType> &inVec );
        // Residual function - hardcoded and AD by hand version
        dType residual_Symbolic( Eigen::MatrixBase<mType> &resVec,
                                  const Eigen::MatrixBase<mType> &solVec);
        dType residual_HandwrittenAdjoint( Eigen::SparseMatrix<dType, Eigen::RowMajor> &Jacobian,
                                 Eigen::MatrixBase<mType> &resVec,
                                 const Eigen::MatrixBase<mType> &solVec);
        template <class vecType>
        dType residual_matFree( const vecType &resVec);
        // Function for loop-wise solver output 
        template <class vecType>
        void writeLoopOutput( const vecType &z, const dType res, const unsigned iteration );
        // Functions to run solver depending on way to determine Jacobian
        void runSolver_WithMatrix( Eigen::MatrixBase<mType> &z, int jacobianOpt, 
                                   dType &res, unsigned &iteration );
        void runSolver_DcoMatrixFree( mType &mz, dType &res, unsigned &iteration );
        
};


template <class mType, class dType> solver<mType, dType>::solver() {
    std::cout << "Construct solver..." << std::endl;
    NminParallel = inputParserObj.getnMinParallel();
    numThreads = inputParserObj.getnumThreads();
    TOL_linsolver = inputParserObj.getTOL_linsolver();
    if ( N >= NminParallel ) {
        Eigen::initParallel();
        Eigen::setNbThreads(numThreads);
    }
}

template <class mType, class dType> solver<mType, dType>::~solver() {
    std::cout << "Close solver..." << std::endl;
}

