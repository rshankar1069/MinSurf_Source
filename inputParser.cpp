/*  Simulation Sciences Laboratory
 *  WS 2019/20
 *  Chenfei Fan, Praveen Mishra, Sankrarasubramanian Ragunathan, Philipp Schleich
 *  Project 1 - "Minimal Surfaces"
 *
 *  Input Parsing class -- source file
 */

#include "inputParser.h"

// Function defined to read the "params.in" file so that the simulation parameters can be read
void input_parser::readFile(std::string filename) {
   std::string line;

   std::ifstream inp_file(filename);
   // std::cout << "++++++++++ Opening the the file \"params.in\" to read the input parameters ++++++++++" << std::endl;
   if(inp_file.is_open()) {
        while(std::getline(inp_file,line)) {
            if (line[0] == '#' || line == "" || line[0] == ' ')
                continue;
            else {
                for(auto i = strtok(&line[0], " "); i != NULL; i = strtok(NULL, " "))
                    tokens.push_back(i);
            }
        }
        inp_file.close();
        // std::cout << "++++++++++ Successfully read the \"params.in\" file ++++++++++" << std::endl;
    }
    else {
        std::cout << "++++++++++ ERROR!! Cannot find \"params.in\" file. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Function defined to set the number of elements along x and y directions of the mesh
void input_parser::setN( ) {
   auto iterator = std::find(tokens.begin(),tokens.end(),"noMeshElements");
   if (iterator != tokens.cend()) {
        int pos = std::distance(tokens.begin(),iterator);
        N = std::stoi(tokens[pos+1]);
    }
    else {
        std::cout << "++++++++++ ERROR!! # of elements along x and y directions not specified. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    } 
}

// Function defined to set the boundary condition at the "bottom" region of the unit square domain
void input_parser::setBCBottom( ) {
    auto iterator = std::find(tokens.begin(),tokens.end(),"BC.BOTTOM");
    if (iterator != tokens.cend()) {
        int pos = std::distance(tokens.begin(),iterator);
        bottom = tokens[pos+1];
    }
    else {
        std::cout << "++++++++++ ERROR!! 'Bottom' boundary condition not specified. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Function defined to set the boundary condition at the "right" region of the unit square domain
void input_parser::setBCRight( ) {
    auto iterator = std::find(tokens.begin(),tokens.end(),"BC.RIGHT");
    if (iterator != tokens.cend()) {
        int pos = std::distance(tokens.begin(),iterator);
        right = tokens[pos+1];
    }
    else {
        std::cout << "++++++++++ ERROR!! 'Right' boundary condition not specified. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Function defined to set the boundary condition at the "top" region of the unit square domain
void input_parser::setBCTop( ) {
    auto iterator = std::find(tokens.begin(),tokens.end(),"BC.TOP");
    if (iterator != tokens.cend()) {
        int pos = std::distance(tokens.begin(),iterator);
        top = tokens[pos+1];
    }
    else {
        std::cout << "++++++++++ ERROR!! 'Top' boundary condition not specified. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Function defined to set the boundary condition at the "left" region of the unit square domain
void input_parser::setBCLeft( ) {
    auto iterator = std::find(tokens.begin(),tokens.end(),"BC.LEFT");
    if (iterator != tokens.cend()) {
        int pos = std::distance(tokens.begin(),iterator);
        left = tokens[pos+1]; 
    }
    else {
        std::cout << "++++++++++ ERROR!! 'Left' boundary condition not specified. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Function defined to set the constants used in the expressions used for the boundary conditions
void input_parser::setConsts( ) {
    auto iterator = std::find(tokens.begin(),tokens.end(),"NUM_CONSTS");
    if (iterator != tokens.cend()) {
        int pos = std::distance(tokens.begin(),iterator);
        int numconsts = std::stoi(tokens[pos+1]);
        if(numconsts != 0) {
            pos = pos+2; // Moving the position to point to the first constant
            for(int i=0; i<2*numconsts; i+=2) {
                consts.insert(std::pair<std::string,float>(tokens[pos+i],std::stof(tokens[pos+i+1])));
            }
        }
        else {
            std::cout << "++++++++++ WARNING!! No external constants specified. Continuing the program" << std::endl;  
        } 
    }
    else {
        std::cout << "++++++++++ WARNING!! No external constants specified. Continuing the program" << std::endl;
    }
}

// Function defined to set the variables used in the expressions used for the boundary conditions
void input_parser::setVars( ) {
    auto iterator = std::find(tokens.begin(),tokens.end(),"NUM_VARS");
    if (iterator != tokens.cend()) {
        int pos = std::distance(tokens.begin(),iterator);
        int numvars = std::stoi(tokens[pos+1]);
        if(numvars != 0) {
            pos = pos+2; // Moving the position to point to the first constant
            for(int i=0; i<2*numvars; i+=2) {
                vars.push_back(tokens[pos+i+1]);
            }
        }
        else {
            std::cout << "++++++++++ ERROR!! Variables not declared in \"params.in\"" << std::endl;  
        } 
    }
    else {
        std::cout << "++++++++++ ERROR!! Variables not declared in \"params.in\"" << std::endl;
    }
}

// Function defined to set the tolerance used for the convergence of the Newton-Raphson Iterations
void input_parser::setTOL( ) {
    auto iterator = std::find(tokens.begin(),tokens.end(),"TOL");
    if (iterator != tokens.cend()) {
        int pos = std::distance(tokens.begin(),iterator);
        TOL = std::stod(tokens[pos+1]);
    }
    else {
        std::cout << "++++++++++ ERROR!! Tolerance for Newton-Raphson Iterations not specified. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Function defined to set the option of using Poisson solution as the initial guess
void input_parser::setPoissonGuess( ) {
    auto iterator = std::find(tokens.begin(),tokens.end(),"usePoissonGuess");
    if (iterator != tokens.cend()) {
        int pos = std::distance(tokens.begin(),iterator);
        usePoissonGuess = std::stoi(tokens[pos+1]);
    }
    else {
        std::cout << "++++++++++ ERROR!! Option to use initial Poisson Solution not specified. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Function defined to set the maximum no. of iterations of the Newton-Raphson Iterations
void input_parser::setmaxIters( ) {
    auto iterator = std::find(tokens.begin(),tokens.end(),"maxIter");
    if (iterator != tokens.cend()) {
        int pos = std::distance(tokens.begin(),iterator);
        maxIter = std::stoi(tokens[pos+1]);
    }
    else {
        std::cout << "++++++++++ ERROR!! Maximum no. of iterations for Newton-Raphson not specified. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Function defined to set the option for using various types of Jacobian for the Newton-Raphson Iterations: hardcoded, ADbyHand, ADbyDCO
void input_parser::setjacOption( ) {
    auto iterator = std::find(tokens.begin(),tokens.end(),"jacobianOpt");
    if (iterator != tokens.cend()) {
        int pos = std::distance(tokens.begin(),iterator);
        jacOption = std::stoi(tokens[pos+1]);
    }
    else {
        std::cout << "++++++++++ ERROR!! Option for the Jacobian not specified. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Function defined to set the option for the frequency at which the output is written
void input_parser::setfileFreq( ) {
    auto iterator = std::find(tokens.begin(),tokens.end(),"fileFreq");
    if (iterator != tokens.cend()) {
        int pos = std::distance(tokens.begin(),iterator);
        fileFreq = std::stoi(tokens[pos+1]);
    }
    else {
        std::cout << "++++++++++ ERROR!! Option for the frequency of file output not specified. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Function defined to return the no. of elements along x and y directions of the mesh
int input_parser::getN( ) {
    return N;
}

// Function defined to return the boundary condition specified at the "bottom" of the unit square domain
std::string input_parser::getBCBottom( ) {
    return bottom;
}

// Function defined to return the boundary condition specified at the "right" of the unit square domain
std::string input_parser::getBCRight( ) {
    return right;
}

// Function defined to return the boundary condition specified at the "top" of the unit square domain
std::string input_parser::getBCTop( ) {
    return top;
}

// Function defined to return the boundary condition specified at the "left" of the unit square domain
std::string input_parser::getBCLeft( ) {
    return left;
}

// Function defined to return the constants used in the expressions for boundary conditions
std::map<std::string,float> input_parser::getConsts( ) {
    return consts;
}

// Function defined to return the variables used in the expressions for boundary conditions
std::vector<std::string> input_parser::getVars( ) {
    return vars;
}

// Function defined to return the tolerance used in the Newton-Raphson Iterations
double input_parser::getTOL( ) {
    return TOL;
}

// Function defined to return the maximum no. of iterations used in the Newton-Raphson Iterations
int input_parser::getmaxIters( ) {
    return maxIter;
}

// Function defined to return the option set for using Poisson solution as the initial guess
int input_parser::getPoissonGuess( ) {
    return usePoissonGuess;
}

// Function defined to return the option set for Jacobian to be used for the Newton-Raphson Iterations
int input_parser::getjacOption( ) {
    return jacOption;
}

// Function defined to return the value of frequency at which the output is written
int input_parser::getfileFreq( ) {
    return fileFreq;
}
