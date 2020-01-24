/*  Simulation Sciences Laboratory
 *  WS 2019/20
 *  Chenfei Fan, Praveen Mishra, Sankarasubramanian Ragunathan, Philipp Schleich
 *  Project 1 - "Minimal Surfaces"
 *
 *  Main execution file
 *  So far: g++ -std=c++14 -O3  minSurf.cpp solver.cpp cartesianGrid.cpp -Wall
 *  @Praveen: could you write a make-file (or CMake, whatever you prefer), that can
 *            -> make (just compile)
 *            -> make run (compile and run the code)
 *            -> make test (compile and run the tests)
 *            -> Tests (you can start with this now I guess)
 *  @someone: -> implement timings (Google benchmark)
 *            -> OpenMP!
 * 
 *  Whatever you do, check out a branch for this "project". I.e. not a personal branch, but a 
 *  branch "inputParser", where you only work on input-stuff or a branch "openMP", where you only 
 *  work on OpenMP
 *  Commit your changes as often as possible inside the branch, and only commit to the master once 
 *  done or once you have achieved a "bigger milestone"
 *  I guess we are doing quite well so far, so we should not have a problem delivering a nice code
 *  at the end!
 */

#include"cartesianGrid.h"
#include"cartesianGrid.cpp"
#include"solver.h"
#include"solver.cpp"

int main(int argc, char *argv[]) {
        
    typedef double dType;
    typedef Eigen::Matrix<dType,-1, 1> mType; // -1 -> Dynamic allocation (for arrays > 16, what 
                                               // we usually have, there is no performance benefit 
                                               // for fixed size allocation
    
    solver<mType, dType> minSurf;
    minSurf.setMesh();
    std::cout << "noGridPoints: " << minSurf.grid.noGridPoints << std::endl;
    if (argc < 2)
        std::cout << "no jacOption" << std::endl;
    else
        minSurf.runSolver(atoi(argv[1]));

    return 0;
}
