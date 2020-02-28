/*  Simulation Sciences Laboratory
 *  WS 2019/20
 *  Chenfei Fan, Praveen Mishra, Sankrarasubramanian Ragunathan, Philipp Schleich
 *  Project 1 - "Minimal Surfaces"
 *
 *  Grid class -- header file
 */

#pragma once


template<class dType, class listType> class cartesianGrid 
{
    public:
        
        // Constructor
        cartesianGrid();
        // Destructor
        ~cartesianGrid();

        // Set information
        void setupGrid( );
       
        // Get information
        dType getCooX(int index);
        dType getCooY(int index);
        
        int noGridPoints;
        dType gridSpacing; // need better name!
        
        listType innerNodeList;
        struct bdryNodes{
            listType bottom, right, top, left;
        };
        bdryNodes bdryNodeList;
        
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // private:                        // Removed Private Dependency for Testing
        
        void setBdryNodes( );
        void setInnerNodes( );
        
        void setNoGridPoints( int N ) { noGridPoints = N; }
        void setGridSpacing( ) { gridSpacing = 1.0 / noGridPoints; } 
                                                                     
        
};


template <class dType, class listType> cartesianGrid<dType, listType>::cartesianGrid() {
    std::cout << "Construct grid..." << std::endl;
}


template <class dType, class listType> cartesianGrid<dType, listType>::~cartesianGrid() {
    std::cout << "Close grid..." << std::endl;
}
