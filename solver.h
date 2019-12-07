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
#include<fstream>
#include"Eigen/Eigen/Sparse"
#include"Eigen/Eigen/Eigenvalues"
#include"Eigen/Eigen/SparseLU"
#include"Eigen/Eigen/SparseCholesky"
#include"Eigen/Eigen/IterativeLinearSolvers"

#include"cartesianGrid.h"
#pragma once

//~#include"inputParser.h" //@Sankar


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
        int N;
        dType h;
        // Boundary conditions
        void applyBC( Eigen::MatrixBase<mType> &inVec );
        
        // Get initial guess by solving -u''=0 in D, u=g on dD via central FD
        void buildPoissonMatrix( Eigen::SparseMatrix<dType> &poissonMatrix );
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
        // Jacobian of minSurf-Operator
        void minSurfJacByHand( Eigen::SparseMatrix<dType> &Jacobian,
                                                            const Eigen::MatrixBase<mType> &inVec);
        //~void minSurfJacByDCO @Chenfei           
        // Residual function
        dType residual( Eigen::MatrixBase<mType> &resVec, const Eigen::MatrixBase<mType> &solVec);
        
};


template <class mType, class dType> solver<mType, dType>::solver() {
    std::cout << "Construct solver..." << std::endl;
}

template <class mType, class dType> solver<mType, dType>::~solver() {
    std::cout << "Close solver..." << std::endl;
}

