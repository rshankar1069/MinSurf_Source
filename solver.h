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
#include<valarray>
#include<cmath>
#include"Eigen/Eigen/Core"
#include"Eigen/Eigen/Sparse"
#include"Eigen/Eigen/IterativeLinearSolvers"
#include"inputParser.h"
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

