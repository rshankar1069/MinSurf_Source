/*  Simulation Sciences Laboratory
 *  WS 2019/20
 *  Chenfei Fan, Praveen Mishra, Sankrarasubramanian Ragunathan, Philipp Schleich
 *  Project 1 - "Minimal Surfaces"
 *
 *  Input Parsing class -- source file
 */

#include "inputParser.h"

input_parser inputParserObj;

// Function defined to read the "params.in" file so that the simulation parameters can be read
void input_parser::readFile(std::string filename)
{
    std::string line;

    std::ifstream inp_file(filename);
    // std::cout << "++++++++++ Opening the the file \"params.in\" to read the input parameters ++++++++++" << std::endl;
    if(inp_file.is_open())
    {
        while(std::getline(inp_file,line))
        {
            if (line[0] == '#' || line == "" || line[0] == ' ')
                continue;
            else
            {
                for(auto i = strtok(&line[0], " "); i != NULL; i = strtok(NULL, " "))
                    tokens.push_back(i);
            }
        }
        inp_file.close();
        // std::cout << "++++++++++ Successfully read the \"params.in\" file ++++++++++" << std::endl;
    }
    else
    {
        std::cout << "++++++++++ ERROR!! Cannot find \"params.in\" file. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Function defined to set the number of elements along x and y directions of the mesh
void input_parser::setN()
{
    auto iterator = std::find(tokens.begin(),tokens.end(),"noMeshElements");
    if (iterator != tokens.cend())
    {
        int pos = std::distance(tokens.begin(),iterator);
        N = std::stoi(tokens[pos+1]);
    }
    else
    {
        std::cout << "++++++++++ ERROR!! # of elements along x and y directions not specified. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    } 
}

// Function defined to set the boundary condition at the "bottom" region of the unit square domain
void input_parser::setBCBottom()
{
    auto iterator = std::find(tokens.begin(),tokens.end(),"BC.BOTTOM");
    if (iterator != tokens.cend())
    {
        int pos = std::distance(tokens.begin(),iterator);
        bottom = tokens[pos+1];
    }
    else
    {
        std::cout << "++++++++++ ERROR!! 'Bottom' boundary condition not specified. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Function defined to set the boundary condition at the "right" region of the unit square domain
void input_parser::setBCRight()
{
    auto iterator = std::find(tokens.begin(),tokens.end(),"BC.RIGHT");
    if (iterator != tokens.cend())
    {
        int pos = std::distance(tokens.begin(),iterator);
        right = tokens[pos+1];
    }
    else
    {
        std::cout << "++++++++++ ERROR!! 'Right' boundary condition not specified. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Function defined to set the boundary condition at the "top" region of the unit square domain
void input_parser::setBCTop()
{
    auto iterator = std::find(tokens.begin(),tokens.end(),"BC.TOP");
    if (iterator != tokens.cend())
    {
        int pos = std::distance(tokens.begin(),iterator);
        top = tokens[pos+1];
    }
    else
    {
        std::cout << "++++++++++ ERROR!! 'Top' boundary condition not specified. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Function defined to set the boundary condition at the "left" region of the unit square domain
void input_parser::setBCLeft()
{
    auto iterator = std::find(tokens.begin(),tokens.end(),"BC.LEFT");
    if (iterator != tokens.cend())
    {
        int pos = std::distance(tokens.begin(),iterator);
        left = tokens[pos+1]; 
    }
    else
    {
        std::cout << "++++++++++ ERROR!! 'Left' boundary condition not specified. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Function defined to set the constants used in the expressions used for the boundary conditions
void input_parser::setConsts()
{
    auto iterator = std::find(tokens.begin(),tokens.end(),"NUM_CONSTS");
    if (iterator != tokens.cend())
    {
        int pos = std::distance(tokens.begin(),iterator);
        int numconsts = std::stoi(tokens[pos+1]);
        if(numconsts != 0)
        {
            pos = pos+2; // Moving the position to point to the first constant
            for(int i=0; i<2*numconsts; i+=2)
            {
                consts.insert(std::pair<std::string,float>(tokens[pos+i],std::stof(tokens[pos+i+1])));
            }
        }
        else
            std::cout << "++++++++++ WARNING!! No external constants specified. Continuing the program" << std::endl;  
    }
    else
    {
        std::cout << "++++++++++ WARNING!! No external constants specified. Continuing the program" << std::endl;
    }
}

// Function defined to set the variables used in the expressions used for the boundary conditions
void input_parser::setVars()
{
    auto iterator = std::find(tokens.begin(),tokens.end(),"NUM_VARS");
    if (iterator != tokens.cend())
    {
        int pos = std::distance(tokens.begin(),iterator);
        int numvars = std::stoi(tokens[pos+1]);
        if(numvars != 0)
        {
            pos = pos+2; // Moving the position to point to the first constant
            for(int i=0; i<2*numvars; i+=2)
            {
                vars.push_back(tokens[pos+i+1]);
            }
        }
        else
            std::cout << "++++++++++ ERROR!! Variables not declared in \"params.in\"" << std::endl;  
    }
    else
    {
        std::cout << "++++++++++ ERROR!! Variables not declared in \"params.in\"" << std::endl;
    }
}

// Function defined to set the tolerance used for the convergence of the Newton-Raphson Iterations
void input_parser::setTOL_Newton()
{
    auto iterator = std::find(tokens.begin(),tokens.end(),"TOL_newton");
    if (iterator != tokens.cend())
    {
        int pos = std::distance(tokens.begin(),iterator);
        TOL_Newton = std::stod(tokens[pos+1]);
    }
    else
    {
        std::cout << "++++++++++ ERROR!! Tolerance for Newton-Raphson iterations not specified. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Function defined to set the tolerance used for the convergence of the respective linear solver 
void input_parser::setTOL_linsolver()
{
    auto iterator = std::find(tokens.begin(),tokens.end(),"TOL_linsolver");
    if (iterator != tokens.cend())
    {
        int pos = std::distance(tokens.begin(),iterator);
        TOL_linsolver = std::stod(tokens[pos+1]);
    }
    else
    {
        std::cout << "++++++++++ ERROR!! Tolerance for linear iterative solver not specified. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Function defined to set the option of using Laplace solution as the initial guess
void input_parser::setInitGuessChoice()
{
    auto iterator = std::find(tokens.begin(),tokens.end(),"initGuessChoice");
    if (iterator != tokens.cend())
    {
        int pos = std::distance(tokens.begin(),iterator);
        initGuessChoice = std::stoi(tokens[pos+1]);
    }
    else
    {
        std::cout << "++++++++++ ERROR!! Option to use initial Laplace Solution not specified. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Function defined to set the maximum no. of iterations of the Newton-Raphson Iterations
void input_parser::setmaxIters()
{
    auto iterator = std::find(tokens.begin(),tokens.end(),"maxIter");
    if (iterator != tokens.cend())
    {
        int pos = std::distance(tokens.begin(),iterator);
        maxIter = std::stoi(tokens[pos+1]);
    }
    else
    {
        std::cout << "++++++++++ ERROR!! Maximum no. of iterations for Newton-Raphson not specified. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Function defined to set the maximum no. of iterations of the Newton-Raphson Iterations
void input_parser::setrelaxNewton()
{
    auto iterator = std::find(tokens.begin(),tokens.end(),"relaxNewton");
    if (iterator != tokens.cend())
    {
        int pos = std::distance(tokens.begin(),iterator);
        relaxNewton = std::stod(tokens[pos+1]);
    }
    else
    {
        std::cout << "++++++++++ ERROR!! Relaxation parameter for Newton iterations not specified. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    }
}


// Function defined to set the option for using various types of Jacobian for the Newton-Raphson Iterations: hardcoded, ADbyHand, ADbyDCO
void input_parser::setjacobianOpt()
{
    auto iterator = std::find(tokens.begin(),tokens.end(),"jacobianOpt");
    if (iterator != tokens.cend())
    {
        int pos = std::distance(tokens.begin(),iterator);
        jacobianOpt = std::stoi(tokens[pos+1]);
    }
    else
    {
        std::cout << "++++++++++ ERROR!! Option for the Jacobian not specified. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Function defined to set the option for the frequency at which the output is written
void input_parser::setfileFreq()
{
    auto iterator = std::find(tokens.begin(),tokens.end(),"fileFreq");
    if (iterator != tokens.cend())
    {
        int pos = std::distance(tokens.begin(),iterator);
        fileFreq = std::stoi(tokens[pos+1]);
    }
    else
    {
        std::cout << "++++++++++ ERROR!! Option for the frequency of file output not specified. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Function to set the option for the no of threads to be used for OpenMP parallelization
void input_parser::setnumThreads()
{
    auto iterator = std::find(tokens.begin(),tokens.end(),"nThreads");
    if (iterator != tokens.cend())
    {
        int pos = std::distance(tokens.begin(),iterator);
        numThreads = std::stoi(tokens[pos+1]);
    }
    else
    {
        std::cout << "++++++++++ ERROR!! Option for the no of threads for parallel run not specified. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Function to set the option for the no of threads to be used for OpenMP parallelization
void input_parser::setnMinParallel()
{
    auto iterator = std::find(tokens.begin(),tokens.end(),"nMinParallel");
    if (iterator != tokens.cend())
    {
        int pos = std::distance(tokens.begin(),iterator);
        nMinParallel = std::stoi(tokens[pos+1]);
    }
    else
    {
        std::cout << "++++++++++ ERROR!! Option for the no of threads for parallel run not specified. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Function to specify the output folder for the vtk files for visualization
void input_parser::setvtkOutLoc()
{
    auto iterator = std::find(tokens.begin(),tokens.end(),"vtkOutLoc");
    if (iterator != tokens.cend())
    {
        int pos = std::distance(tokens.begin(),iterator);
        vtkOutLoc = tokens[pos+1];
    }
    else
    {
        std::cout << "++++++++++ ERROR!! Option for the vtk file output location not specified. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Function to specify the output folder for the residual files for visualization
void input_parser::setresOutLoc()
{
    auto iterator = std::find(tokens.begin(),tokens.end(),"resOutLoc");
    if (iterator != tokens.cend())
    {
        int pos = std::distance(tokens.begin(),iterator);
        resOutLoc = tokens[pos+1];
    }
    else
    {
        std::cout << "++++++++++ ERROR!! Option for the residual file output location not specified. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Function defined to return the no. of elements along x and y directions of the mesh
int input_parser::getN() {
    return N;
}

// Function defined to return the boundary condition specified at the "bottom" of the unit square domain
std::string input_parser::getBCBottom() {
    return bottom;
}

// Function defined to return the boundary condition specified at the "right" of the unit square domain
std::string input_parser::getBCRight() {
    return right;
}

// Function defined to return the boundary condition specified at the "top" of the unit square domain
std::string input_parser::getBCTop() {
    return top;
}

// Function defined to return the boundary condition specified at the "left" of the unit square domain
std::string input_parser::getBCLeft() {
    return left;
}

// Function defined to return the constants used in the expressions for boundary conditions
std::map<std::string,float> input_parser::getConsts() {
    return consts;
}

// Function defined to return the variables used in the expressions for boundary conditions
std::vector<std::string> input_parser::getVars() {
    return vars;
}

// Function defined to return the tolerance used in the Newton-Raphson Iterations
double input_parser::getTOL_Newton() {
    return TOL_Newton;
}

// Function defined to return the tolerance used in the linear iterative solver
double input_parser::getTOL_linsolver() {
    return TOL_linsolver;
}

// Function defined to return the maximum no. of iterations used in the Newton-Raphson Iterations
int input_parser::getmaxIters() {
    return maxIter;
}

// Function defined to return the relaxation parameter for the Newton procedure 
double input_parser::getrelaxNewton() {
    return relaxNewton;
}

// Function defined to return the option set for using Laplace solution as the initial guess
int input_parser::getInitGuessChoice() {
    return initGuessChoice;
}

// Function defined to return the option set for Jacobian to be used for the Newton-Raphson Iterations
int input_parser::getjacobianOpt() {
    return jacobianOpt;
}

// Function defined to return the value of frequency at which the output is written
int input_parser::getfileFreq() {
    return fileFreq;
}

// Function defined to return the no of threads to be used for parallel execution of the program
int input_parser::getnumThreads() {
    return numThreads;
}

// Function defined to return the no of threads to be used for parallel execution of the program
int input_parser::getnMinParallel() {
    return nMinParallel;
}

// Function defined to return the vtk file output location
std::string input_parser::getvtkOutLoc() {
    return vtkOutLoc;
}

// Function defined to return the residual file output location
std::string input_parser::getresOutLoc() {
    return resOutLoc;
}
