#include "inputParser.h"

void input_parser::readFile(std::string filename)
{
    std::string line;

    std::ifstream inp_file(filename);
    std::cout << "++++++++++ Opening the the file \"params.in\" to read the input parameters ++++++++++" << std::endl;
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
        std::cout << "++++++++++ Successfully read the \"params.in\" file ++++++++++" << std::endl;
    }
    else
    {
        std::cout << "++++++++++ ERROR!! Cannot find \"params.in\" file. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    }
}

void input_parser::setN()
{
    auto iterator = std::find(tokens.begin(),tokens.end(),"N");
    if (iterator != tokens.cend())
    {
        int pos = std::distance(tokens.begin(),iterator);
        N = std::stoi(tokens[pos+1]);
        std::cout << "N\t==\t" << N << std::endl; 
    }
    else
    {
        std::cout << "++++++++++ ERROR!! # of elements along x and y directions not specified. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    } 
}

void input_parser::setBCBottom()
{
    auto iterator = std::find(tokens.begin(),tokens.end(),"BC.BOTTOM");
    if (iterator != tokens.cend())
    {
        int pos = std::distance(tokens.begin(),iterator);
        bottom = tokens[pos+1];
        std::cout << "BC.BOTTOM\t==\t" << bottom << std::endl; 
    }
    else
    {
        std::cout << "++++++++++ ERROR!! 'Bottom' boundary condition not specified. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    }
}

void input_parser::setBCRight()
{
    auto iterator = std::find(tokens.begin(),tokens.end(),"BC.RIGHT");
    if (iterator != tokens.cend())
    {
        int pos = std::distance(tokens.begin(),iterator);
        right = tokens[pos+1];
        std::cout << "BC.RIGHT\t==\t" << right << std::endl; 
    }
    else
    {
        std::cout << "++++++++++ ERROR!! 'Right' boundary condition not specified. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    }
}

void input_parser::setBCTop()
{
    auto iterator = std::find(tokens.begin(),tokens.end(),"BC.TOP");
    if (iterator != tokens.cend())
    {
        int pos = std::distance(tokens.begin(),iterator);
        top = tokens[pos+1];
        std::cout << "BC.TOP\t==\t" << top << std::endl; 
    }
    else
    {
        std::cout << "++++++++++ ERROR!! 'Top' boundary condition not specified. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    }
}

void input_parser::setBCLeft()
{
    auto iterator = std::find(tokens.begin(),tokens.end(),"BC.LEFT");
    if (iterator != tokens.cend())
    {
        int pos = std::distance(tokens.begin(),iterator);
        left = tokens[pos+1];
        std::cout << "BC.LEFT\t==\t" << left << std::endl; 
    }
    else
    {
        std::cout << "++++++++++ ERROR!! 'Left' boundary condition not specified. Exiting the program ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    }
}

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

int input_parser::getN()
{
    return N;
}

std::string input_parser::getBCBottom()
{
    return bottom;
}

std::string input_parser::getBCRight()
{
    return right;
}

std::string input_parser::getBCTop()
{
    return top;
}

std::string input_parser::getBCLeft()
{
    return left;
}

std::map<std::string,float> input_parser::getConsts()
{
    return consts;
}

std::vector<std::string> input_parser::getVars()
{
    return vars;
}
