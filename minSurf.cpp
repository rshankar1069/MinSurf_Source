/*  Simulation Sciences Laboratory
 *  WS 2019/20
 *  Chenfei Fan, Praveen Mishra, Sankrarasubramanian Ragunathan, Philipp Schleich
 *  Project 1 - "Minimal Surfaces"
 *
 *  Main execution file
*/

#include"inputParser.h"
#include"inputParser.cpp"
#include"cartesianGrid.h"
#include"cartesianGrid.cpp"
#include"solver.h"
#include"solver.cpp"
#include"postProcessor.h"
#include"postProcessor.cpp"


int main() {

    // Change dType if you want to run the simulation in float or double
    typedef double dType;
    typedef Eigen::Matrix<dType,-1, 1> mType; // -1 -> Dynamic allocation (for arrays > 16, what 
                                               // we usually have, there is no performance benefit 
                                               // for fixed size allocation
   
    // Create instance of solver 
    solver<mType, dType> minSurf;
    // Setup the mesh
    minSurf.setMesh();
    // Output number of grid points
    std::cout << "-------------------------------------------------------------" << std::endl;
    std::cout << "N: " << minSurf.grid.noGridPoints << "\t number of gridpoints: "
              << minSurf.grid.noGridPoints*minSurf.grid.noGridPoints << std::endl;
    std::cout << "-------------------------------------------------------------" << std::endl;
    // Run the solver, including pre- & postprocessing according to input file
    minSurf.runSolver();

}
