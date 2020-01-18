/*  Simulation Sciences Laboratory
 *  WS 2019/20
 *  Chenfei Fan, Praveen Mishra, Sankarasubramanian Ragunathan, Philipp Schleich
 *  Project 1 - "Minimal Surfaces"
 *
 *  Grid class -- header file
 */

#pragma once

// Here we want some further explanations...

template<class dType, class listType> class cartesianGrid 
{
    public:
        
        // Constructor
        cartesianGrid();
        // Destructor
        ~cartesianGrid();

        // Set information
        void setupGrid( );
                  // @Sankar: you could define this function for the marriage 
                  // with the input parser! (then inside, call setNoGridPoints etc. to
                  // set everything important up)
        
        // Get information
        dType getCooX(int index);
        dType getCooY(int index);
        
        // Grid parameters
        int noGridPoints;
        dType gridSpacing; // need better name!
        
        // Node lists containing inner and boundary nodes
        listType innerNodeList;
        struct bdryNodes{
            listType bottom, right, top, left;
        };
        bdryNodes bdryNodeList;
        
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    private:
        
        // Methods to set up node lists for structured grid
        void setBdryNodes( );
        void setInnerNodes( );
        
        // Methods to set up grid parameters
        void setNoGridPoints( int N ) { noGridPoints = N; }
        void setGridSpacing( ) { gridSpacing = 1.0 / noGridPoints; } // make sure to only use when 
                                                                     // noGridPoints is set!
        
};


template <class dType, class listType> cartesianGrid<dType, listType>::cartesianGrid() {
    std::cout << "Construct grid..." << std::endl;
}


template <class dType, class listType> cartesianGrid<dType, listType>::~cartesianGrid() {
    std::cout << "Close grid..." << std::endl;
}
