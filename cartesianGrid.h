/*  Simulation Sciences Laboratory
 *  WS 2019/20
 *  Chenfei Fan, Praveen Mishra, Sankrarasubramanian Ragunathan, Philipp Schleich
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
        
        // We could also define getters for the next two variables, but I don't know if that won't 
        // be too inconvenient later
        int noGridPoints;
        dType gridSpacing; // need better name!
        
        listType innerNodeList;
        struct bdryNodes{
            listType bottom, right, top, left;
        };
        bdryNodes bdryNodeList;
        
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    private:
        
        void setBdryNodes( );
        void setInnerNodes( );
        
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
