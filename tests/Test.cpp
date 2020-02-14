/*  Simulation Sciences Laboratory
 *  WS 2019/20
 *  Chenfei Fan, Praveen Mishra, Sankarasubramanian Ragunathan, Philipp Schleich
 *  Project 1 - "Minimal Surfaces"
 *
 *  TEST FIXTURES
*/ 
 
#include "../inputParser.h"  
#include "../inputParser.cpp" 
#include "../solver.h"  
#include "../solver.cpp"  
#include "../cartesianGrid.h"  
#include "../cartesianGrid.cpp"
#include "gtest/gtest.h"


input_parser SineTest("../testcases/sine/params.in");
input_parser inputParserobj("../testcases/sine/params.in");
//input_parser ScherkTest("../testcase/scherk/params.in");


// TEST Initial Predefined Conditions For Sine

TEST(Testinputs, InputParsalSine){
	
	EXPECT_TRUE(SineTest.getN()==300);       // Number of Mesh Elements= 300
	EXPECT_GT(SineTest.getmaxIters(), 100);
	EXPECT_EQ(SineTest.getnumThreads(), 4);
    EXPECT_LT(100, SineTest.getnMinParallel());
	EXPECT_NE(1e-5, SineTest.getTOL_Newton());    //  Maximumum Tolerance - 1e-5
	EXPECT_EQ( SineTest.getrelaxNewton() , 0.9);
	
	ASSERT_TRUE(SineTest.getPoissonGuess()== 1);  // Fatal Test. Default solver is by using initial Poisson Guess(1)
	ASSERT_FALSE(SineTest.getTOL_linsolver() ==1e-5);
	
	
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

// Test to check for Cartesian Grid Formation. 
//Calling all the functions within class Cartesian Grid. Manually Providing value for 'N' and noGridPoints

TEST(Testgrid, BoundaryNodes){
	
int bottom_most,left_most,top_most,right_most, i=0;		
cartesianGrid<double, std::vector<int>> gridT;
int N= SineTest.getN();
int noGridPoints=N;
gridT.setNoGridPoints(N);
gridT.setGridSpacing(); 
gridT.setBdryNodes();
gridT.setInnerNodes();
std::cout << N << std::endl;
	int bot,rig,to,lef;
    for(auto& i: gridT.bdryNodeList.bottom)
    bottom_most=i; 
    for(auto& i: gridT.bdryNodeList.right)
    right_most=i;
    for(auto& i: gridT.bdryNodeList.top)
    top_most=i;
    for(auto& i: gridT.bdryNodeList.left)
    left_most=i;

int size= gridT.innerNodeList.size();

EXPECT_TRUE(size==100);
EXPECT_LT(left_most, top_most);
EXPECT_GT(top_most, right_most);
EXPECT_GT(right_most, bottom_most);
EXPECT_LT(bottom_most,top_most);
EXPECT_LT(bottom_most, left_most); 


EXPECT_EQ(gridT.noGridPoints, 300);
}



TEST(TestSolver, BC){
	
typedef Eigen::Matrix<double,-1, 1> mType;

solver<mType, double>solution;
EXPECT_EQ(1 ,2);
mType invec;
invec << 1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0; 
std::cout << invec[4] << std::endl;
solution.applyBC(invec);

std::cout << invec[4] << std::endl;



}




int main(int argc, char* argv[]) {
	
  ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
	


}

