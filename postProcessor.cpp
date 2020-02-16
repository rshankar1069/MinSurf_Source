/*  Simulation Sciences Laboratory
 *  WS 2019/20
 *  Chenfei Fan, Praveen Mishra, Sankrarasubramanian Ragunathan, Philipp Schleich
 *  Project 1 - "Minimal Surfaces"
 *
 *  Post-processing class
 */

#include "postProcessor.h"

// Function to write the structured VTK file of the final solution for visualization
template<class mType,class dType>
void structuredGridWriter(int iterationIndex, mType z)
{
    // Input parser object to get the number of nodes
    int N = inputParserObj.getN();

    // Remove directory for every new simulation
    std::stringstream removedir;
    removedir << "rm -R " << inputParserObj.getvtkOutLoc();
    struct stat info;
    if(iterationIndex == inputParserObj.getfileFreq() || iterationIndex == 1 && stat(removedir.str().c_str(),&info) != 0){
        system(removedir.str().c_str());
    }

    // Writing the solution output into .vts file
    std::stringstream ss,makedir;
    makedir << "mkdir -p " << inputParserObj.getvtkOutLoc();
    system(makedir.str().c_str());
    ss << inputParserObj.getvtkOutLoc() << "vtk" << iterationIndex << ".dat";
    std::string filename = ss.str();

    std::ofstream file;
    file.open(filename,std::ios::out | std::ios::trunc);
    file.precision(15);
    for (auto i=0; i<z.size(); i++) file << z[i] << std::endl;
    file.close();

    // Reading the output file to process the data into .vtk format
    std::ifstream vtkf;
    vtkf.open(filename);

    dType h = 1.0/N;

    if(!vtkf.is_open())
    {
        std::cout << "++++++++++ Error in opening the file !! Failed to post-process the data ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    }
    else
    {
        // Storing the output data in a vector
        std::vector<dType> data;
        std::string str;
        while (std::getline(vtkf,str))
        {
            if(str.size()>0)
                data.push_back(std::stod(str));
        }

        // Create a grid
        vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
        
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        
        unsigned k=0;
        for(unsigned int j=0;j<N;j++)
            for(unsigned int i=0;i<N;i++)
            {
                points->InsertNextPoint(i*h,j*h,data[k]);
                k++;
            }
        
        // Specify the dimensions of the grid
        structuredGrid->SetDimensions(N,N,1);
        structuredGrid->SetPoints(points);
        
        // Write file
        vtkSmartPointer<vtkXMLStructuredGridWriter> writer =
            vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
        
        std::stringstream vtkfile;
        vtkfile << inputParserObj.getvtkOutLoc() << "vtk_" << iterationIndex << ".vts";
        writer->SetFileName(vtkfile.str().c_str());
        writer->SetInputData(structuredGrid);
        writer->Write();
    }
}

// Function defined to write the residual values into a csv file
template <class dType>
void residualWriter(int iterationIndex, dType res)
{
    int N = inputParserObj.getN();

    // Remove directory for every new simulation
    std::stringstream removedir;
    removedir << "rm -R " << inputParserObj.getresOutLoc();
    struct stat info;
    if(iterationIndex == inputParserObj.getfileFreq() || iterationIndex == 1 && stat(removedir.str().c_str(),&info) != 0){
        system(removedir.str().c_str());
    }

    // Writing the residual to a csv file
    std::stringstream ss,makedir;
    makedir << "mkdir -p " << inputParserObj.getresOutLoc();
    system(makedir.str().c_str());
    
    ss << inputParserObj.getresOutLoc() << "residual" << ".dat";
    std::string filename = ss.str();

    std::ofstream resfile;
    resfile.open(filename,std::ios::out | std::ios::app);
    if(iterationIndex == 1)
    {
        resfile << "\"Iteration Index\"" << "," << "\"Residual\"" << std::endl;
        resfile << iterationIndex << "," << res << std::endl;
    }
    else
    {
        resfile << iterationIndex << "," << res << std::endl;
    }
    resfile.close();
}
