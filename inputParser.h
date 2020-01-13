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
        double TOL;                                         // Tolerance value for the Newton-Raphson Iterations
        int usePoissonGuess;                                // Option to choose Poisson Guess for the initial condition
        int maxIter;                                        // Maximum number of iterations
        int jacobianOpt;                                    // Option to choose the Jacobian to be used in solving the problem
        int fileFreq;                                       // Option to choose the frequency in which the files are written

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
        void setTOL();
        void setPoissonGuess();
        void setjacobianOpt();
        void setmaxIters();
        void setfileFreq();
        
        // Getter Function Definitions
        int getN();
        std::string getBCBottom();
        std::string getBCRight();
        std::string getBCTop();
        std::string getBCLeft();
        std::map<std::string,float> getConsts();
        std::vector<std::string> getVars();
        double getTOL();
        int getPoissonGuess();
        int getjacobianOpt();
        int getmaxIters();
        int getfileFreq();

        // Constructer definition to call Setter Functions
        input_parser(std::string filename="../params.in")
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
            setTOL();
            setPoissonGuess();
            setmaxIters();
            setjacobianOpt();
            setfileFreq();
        }
};
#endif //INPUTPARSER_H_INCLUDED
