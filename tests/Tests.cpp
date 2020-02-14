#include "../inputParser.h"  
#include "../inputParser.cpp" 
#include "../solver.h"  
#include "../solver.cpp"  
#include "../cartesianGrid.h"  
#include "../cartesianGrid.cpp"
#include "gtest/gtest.h"


input_parser SineTest("../testcases/sine/params.in");
//input_parser ScherkTest("../testcase/scherk/params.in");


// TEST Initial Predefined Conditions For Sine

TEST(Testinputs, InputParsalSine){
	
	EXPECT_TRUE(SineTest.getN()==300);       // Number of Mesh Elements= 300
	EXPECT_GT(SineTest.getmaxIters(), 100);
	EXPECT_EQ(SineTest.getnumThreads(), 4);
    EXPECT_LT(100, SineTest.getnMinParallel());
	EXPECT_NE(1e-5, SineTest.getTOL_Newton());    //  Maximumum Tolerance - 1e-5
	EXPECT_EQ( SineTest.getrelaxNewton() , 0.9);
	
	//ASSERT_TRUE(SineTest.getPoissonGuess(), 1);  // Fatal Test. Default solver is by using initial Poisson Guess(1)
	//ASSERT_FALSE(SineTest.getTOL_linsolver() , 1e-5);
	
	
}

/*
TEST(Testinputs, InputParsalScherk){
	
	EXPECT_TRUE(ScherkTest.getN()==50);      
	EXPECT_GT(ScherkTest.getmaxIters(), 100);
	EXPECT_EQ(ScherkTest.getnumThreads(), 4);
    EXPECT_LT(100, ScherkTest.getnMinParallel());
	EXPECT_NE(1e-6, ScherkTest.getTOL_Newton());    //  Maximumum Tolerance - 1e-6
	EXPECT_EQ( ScherkTest.getrelaxNewton() , 0.9);
	
	//ASSERT_TRUE(ScherkTest.getPoissonGuess(), 1);  // Fatal Test. Default solver is by using initial Poisson Guess(1)
	//ASSERT_FALSE(ScherkTest.getTOL_linsolver() , 1e-5);
	
	
}


*/


TEST(Testgrid, BoundaryNodes){
		
cartesianGrid<int, std::list<int>> gridT;
int N= SineTest.getN();
std::cout << N << std::endl;
/*	
   
	int bot,rig,to,lef;
    for(auto& i: gridT.bdryNodeList.bottom)
    bottom_most=i; 
    for(auto& i: gridT.bdryNodeList.right)
    right_most=i;
    for(auto& i: gridT.bdryNodeList.top)
    top_most=i;
    for(auto& i: gridT.bdryNodeList.left)
    left_most=i;

EXPECT_GT(left_most, top_most);
EXPECT_GT(top_most, right_most);
EXPECT_GT(right_most, bottom_most);
EXPECT_LT(bottom_most,top_most);
EXPECT_LT(bottom_most, left_most); 

*/

EXPECT_GT(gridT.noGridPoints, 300);
}


TEST(Testgrid, Innernodelist){
	
	
	cartesianGrid<int, std::list<int>> testgrid;
    int k = testgrid.innerNodeList.size();
	
	EXPECT_EQ(k, 100);
        
}



// Debug TestCase
TEST(TestCase1, TestName1) {
  EXPECT_EQ(120, 120);
  EXPECT_TRUE(true);
}


TEST(Solver, Setmesh){	
solver<int, int> Solution;
//Solution.buildLaplaceMatrix( Eigen::SparseMatrix<dType, Eigen::RowMajor> &laplaceMatrix );
//Solution.laplaceMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
	
/*
build friend class
-> pass by reference -> use extern input parser	
	
	
}
*/
int main(int argc, char* argv[]) {
	
  ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
	


}

