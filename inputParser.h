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
        
        ATMSP<float> parser;                                // Expression parser for parsing the boundary conditions
        ATMSB<float> byteCode;                              // Expression parser respresented in terms of byte code
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
        
        // Getter Function Definitions
        int getN();
        std::string getBCBottom();
        std::string getBCRight();
        std::string getBCTop();
        std::string getBCLeft();
        std::map<std::string,float> getConsts();
        std::vector<std::string> getVars();

        // Constructer definition to call Setter Functions
        input_parser(std::string filename="params.in")
        {
            // Read the file initially
            #pragma once
            readFile(filename);
            setN();
            setBCBottom();
            setBCRight();
            setBCTop();
            setBCLeft();
            setConsts();
            setVars();
        }
};
#endif //INPUTPARSER_H_INCLUDED
