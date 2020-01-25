/*  Simulation Sciences Laboratory
 *  WS 2019/20
 *  Chenfei Fan, Praveen Mishra, Sankarasubramanian Ragunathan, Philipp Schleich
 *  Project 1 - "Minimal Surfaces"
 *
 *  Solver class -- header file
 */

#include<iostream>
#include<stdio.h>
#include<vector>
#include<valarray>
#include<cmath>
#include"Eigen/Eigen/Sparse"
#include"Eigen/Eigen/Eigenvalues"
#include"Eigen/Eigen/SparseLU"
#include"Eigen/Eigen/SparseCholesky"
#include"Eigen/Eigen/IterativeLinearSolvers"

#include"cartesianGrid.h"
#pragma once


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
        
        void setMesh( );
        void runSolver( );
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    private:
        // Private variables
        int N; // number of gridpoints
        dType h; // grid spacing
        int jacOption; // switch for options to determine Jacobian
        int numThreads = 16; // number of threads for OpenMP multithreading
        int NminParallel = 10;
 
        // Private methods
        
        // Boundary conditions
        void applyBC( Eigen::MatrixBase<mType> &inVec );
        
        // Get initial guess by solving -u''=0 in D, u=g on dD via central FD
        void buildPoissonMatrix( Eigen::SparseMatrix<dType, Eigen::RowMajor> &poissonMatrix );
        void getInitGuess( Eigen::MatrixBase<mType> &zE );
        
        // Helper functions for FD stencils
        inline dType getDx( const Eigen::MatrixBase<mType> &inVec, const int index );
        inline dType getDy( const Eigen::MatrixBase<mType> &inVec, const int index );
        inline dType getDxx( const Eigen::MatrixBase<mType> &inVec, const int index );
        inline dType getDyy( const Eigen::MatrixBase<mType> &inVec, const int index );
        inline dType getDxy( const Eigen::MatrixBase<mType> &inVec, const int index );
        
        // minSurf-Operator
        void minSurfOperator( Eigen::MatrixBase<mType> &outVec,
                        const Eigen::MatrixBase<mType> &inVec );
        // Jacobian of minSurf-Operator - hardcoded and AD by hand version
        void minSurfJac_HardCoded( Eigen::SparseMatrix<dType, Eigen::RowMajor> &Jacobian,
                               const Eigen::MatrixBase<mType> &inVec);
        void minSurfJac_ADByHand( Eigen::SparseMatrix<dType, Eigen::RowMajor> &Jacobian,
                                 Eigen::MatrixBase<mType> &outVec,
                                 const Eigen::MatrixBase<mType> &inVec );
        // Residual function - hardcoded and AD by hand version
        dType residual_HardCoded( Eigen::MatrixBase<mType> &resVec,
                                  const Eigen::MatrixBase<mType> &solVec);
        dType residual_ADByHand( Eigen::SparseMatrix<dType, Eigen::RowMajor> &Jacobian, 
                                 Eigen::MatrixBase<mType> &resVec,
                                 const Eigen::MatrixBase<mType> &solVec);
                                 
        // Functions to run solver depending on way to determine Jacobian
        void runSolver_HardCoded( );
        void runSolver_ADByHand( );
        
};


template <class mType, class dType> solver<mType, dType>::solver() {
    std::cout << "Construct solver..." << std::endl;
    if ( N >= NminParallel ) {
        Eigen::initParallel();
        Eigen::setNbThreads(numThreads);
    }
}

template <class mType, class dType> solver<mType, dType>::~solver() {
    std::cout << "Close solver..." << std::endl;
}

