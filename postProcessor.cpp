/*  Simulation Sciences Laboratory
 *  WS 2019/20
 *  Chenfei Fan, Praveen Mishra, Sankrarasubramanian Ragunathan, Philipp Schleich
 *  Project 1 - "Minimal Surfaces"
 *
 *  Post-processing class
 */

#include "postProcessor.h"

// Function to write the structured VTK file of the final solution for visualization
template<class dType>
void structuredGridWriter(std::string filename,int iterationIndex) {
    // Reading the output file to process the data into .vtk format
    std::ifstream file;
    file.open(filename);
    // Input parser object to get the number of nodes
    input_parser inputParserObj;
    int N = inputParserObj.getN();

    dType h = 1.0/N;

    if(!file.is_open()) {
        std::cout << "++++++++++ Error in opening the file !! Failed to post-process the data ++++++++++" << std::endl;
        exit(EXIT_FAILURE);
    }
    else {
        // Storing the output data in a vector
        std::vector<dType> data;
        std::string str;
        while (std::getline(file,str)) {
            if(str.size()>0)
                data.push_back(std::stod(str));
        }

        // Create a grid
        vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
        
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        
        unsigned k=0;
        for(unsigned int j=0;j<N;j++) {
            for(unsigned int i=0;i<N;i++) {
                points->InsertNextPoint(i*h,j*h,data[k]);
                k++;
            }
        }
        
        // Specify the dimensions of the grid
        structuredGrid->SetDimensions(N,N,1);
        structuredGrid->SetPoints(points);
        
        // Write file
        vtkSmartPointer<vtkXMLStructuredGridWriter> writer =
            vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
        
        char vtsfile[50];
        sprintf(vtsfile,"../Output/Solution/output_%d.vts",iterationIndex);
        writer->SetFileName(vtsfile);
        writer->SetInputData(structuredGrid);
        writer->Write();
    }
}
