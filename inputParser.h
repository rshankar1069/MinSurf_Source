/*  Simulation Sciences Laboratory
 *  WS 2019/20
 *  Chenfei Fan, Praveen Mishra, Sankrarasubramanian Ragunathan, Philipp Schleich
 *  Project 1 - "Minimal Surfaces"
 *
 *  Input Parsing class -- header file
 */

#ifndef INPUTPARSER_H_INCLUDED
#define INPUTPARSER_H_INCLUDED

/* Required include files for the input file parsing */
#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <map>
#include <vector>
#include "atmsp.h"

class input_parser
{
    private:
        int N;                                              // # of elements along x and y
        std::string bottom,right,top,left;                  // BC expressions along the edges of the domain
        std::map<std::string,float> consts;                 // Extra constants used in defining the expressions
        std::vector<std::string> vars;                      // Variables defined in the expressions
        std::vector<std::string> tokens;
        std::string varnames;                               // Temporary variable to store variable names
        double TOL_Newton;                                  // Tolerance value for the Newton-Raphson Iterations
        double TOL_linsolver;                               // Tolerance value for linear solver
        int initGuessChoice;                                // Option to choose initial condition
        int maxIter;                                        // Maximum number of iterations
        double relaxNewton;                                 // relaxation parameter for Newton procedure
        int jacobianOpt;                                    // Option to choose the Jacobian to be used in solving the problem
        int fileFreq;                                       // Option to choose the frequency in which the files are written
        int numThreads;                                     // Option to choose the no of threads to use for OpenMP
        int nMinParallel;                                   // Option to choose minimal number of gridpoints to execute in parallel
        int scherkCheck;                                    // Option to compare the FD solution with Analytical Scherk equation
        std::string vtkOutLoc;                              // Output Location to store the vtk files for postprocessing
        std::string resOutLoc;                              // Output Location to store the residual files for postprocessing

    public:
        // Function to read the file and initialize the mesh parameters
        void readFile(std::string filename);
        
        // Setter Function Definitions
        void setN();
        void setBCBottom();
        void setBCRight();
        void setBCTop();
        void setBCLeft();
        void setConsts();
        void setVars();
        void setTOL_Newton();
        void setTOL_linsolver();
        void setInitGuessChoice();
        void setjacobianOpt();
        void setmaxIters();
        void setrelaxNewton();
        void setfileFreq();
        void setnumThreads();
        void setnMinParallel();
        void setscherkCheck();
        void setvtkOutLoc();
        void setresOutLoc();
        
        // Getter Function Definitions
        int getN();
        std::string getBCBottom();
        std::string getBCRight();
        std::string getBCTop();
        std::string getBCLeft();
        std::map<std::string,float> getConsts();
        std::vector<std::string> getVars();
        double getTOL_Newton();
        double getTOL_linsolver();
        int getInitGuessChoice();
        int getjacobianOpt();
        int getmaxIters();
        double getrelaxNewton();
        int getfileFreq();
        int getnumThreads();
        int getnMinParallel();
        int getscherkCheck();
        std::string getvtkOutLoc();
        std::string getresOutLoc();

        // Constructer definition to call Setter Functions
        input_parser(std::string filename="params.in")
        {
            // Read the file initially
            readFile(filename);
            setN();
            setBCBottom();
            setBCRight();
            setBCTop();
            setBCLeft();
            setConsts();
            setVars();
            setTOL_Newton();
            setTOL_linsolver();
            setInitGuessChoice();
            setmaxIters();
            setrelaxNewton();
            setjacobianOpt();
            setfileFreq();
            setnumThreads();
            setnMinParallel();
            setscherkCheck();
            setvtkOutLoc();
            setresOutLoc();    
        }
};
#endif //INPUTPARSER_H_INCLUDED
