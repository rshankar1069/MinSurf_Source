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

typedef Eigen::Matrix<double,-1, 1> mType;

//#######################################################################################

// TEST input Parser for Sine Testcase

TEST(Testinputs, InputParserSine){
	
	EXPECT_TRUE(inputParserObj.getN()==300);       // Number of Mesh Elements= 300
	EXPECT_GT(inputParserObj.getmaxIters(), 100);
	EXPECT_EQ(inputParserObj.getnumThreads(), 4);
    EXPECT_LT(100, inputParserObj.getnMinParallel());
	EXPECT_NE(1e-5, inputParserObj.getTOL_Newton());    //  Maximumum Tolerance - 1e-5
	EXPECT_EQ( inputParserObj.getrelaxNewton() , 0.9);
	
//	ASSERT_TRUE(inputParserObj.getPoissonGuess()== 1);  // Fatal Test. Default solver is by using initial Poisson Guess(1)
	ASSERT_FALSE(inputParserObj.getTOL_linsolver() ==1e-5);
	
	
}
/*
TEST(Testinputs, InputParserScherk){
	
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

//##########################################################################################################
// Test to check for Cartesian Grid Formation. 
//Calling all the functions within class Cartesian Grid. Manually Providing value for 'N' and noGridPoints

TEST(Testgrid, BoundaryNodes){
	
	int bottom_most,left_most,top_most,right_most, i=0;		
	cartesianGrid<double, std::vector<int>> gridT;
	int N= inputParserObj.getN();
	int noGridPoints=N;
	gridT.setNoGridPoints(N);
	gridT.setGridSpacing(); 
	gridT.setBdryNodes();
	gridT.setInnerNodes();
	std::cout << N << std::endl;

    for(auto& i: gridT.bdryNodeList.bottom)
    bottom_most=i; 
    for(auto& i: gridT.bdryNodeList.right)
    right_most=i;
    for(auto& i: gridT.bdryNodeList.top)
    top_most=i;
    for(auto& i: gridT.bdryNodeList.left)
    left_most=i;


}
// Testing for Lexicographical ordering

TEST(Testgrid, InnerNodes){
	
int end;	
	
	cartesianGrid<double, std::vector<int>> gridT2;  // Test instance of class cartesianGrid
	int N= inputParserObj.getN();
	int noGridPoints=N;
	gridT2.setNoGridPoints(N);                       // setting up grid   
	gridT2.setGridSpacing(); 
	gridT2.setBdryNodes();                           // sorting
	gridT2.setInnerNodes();

for(auto& it: gridT2.innerNodeList)
	it=end; 		
int size= gridT2.innerNodeList.size();
EXPECT_EQ(gridT2.noGridPoints, 300);
EXPECT_GT(size , 100);
EXPECT_EQ(end, 0);            
}
// ##########################################################################################
// Testing for boundary conditions with zero matrix

TEST(TestSolver, BoundaryConditions){
	int f1,f2;	
	int N=inputParserObj.getN();
	solver<mType, double>solution;                      // Test instance of class solver
	mType testvec = mType::Zero(N*N);                  
	solution.setMesh();
	solution.applyBC(testvec);                          // Applying BC to test vector
	for(int i=0; i<N; i++){
	if (i ==0)
	f1= testvec[i];
else
	if (i==N-1)
		f2= testvec[i];
}	

	EXPECT_EQ(f1,1);                               
	EXPECT_EQ(f2, 1);

}


TEST (TestSolver , LaplaceMatrix ){
	int N=3;	
	solver<mType, double>solution2;	
	Eigen::SparseMatrix<double, Eigen::RowMajor> testmatrix(N*N, N*N);
	solution2.buildLaplaceMatrix(testmatrix);
	
	EXPECT_EQ(testmatrix.rows(), 10);	

}
/*
TEST (TestSolver,residual_HandwrittenAdjoint){
	
	Eigen::SparseMatrix<dType, Eigen::RowMajor> Jacobian(N*N, N*N);
	mType test = mType::Zero(N*N);
	const mTypey reference= mType::Zero(N*N);
	solver<mType,double>solution3;
	double test_residual=solution3.residual_HandwrittenAdjoint(Jacobian,test,reference);
	
	EXPECT_TRUE( test_residual == 1);
		
		

}
*/

int main(int argc, char* argv[]) {
	// std::cout << inputParserObj.getN() << std::endl; // Debugging Input Parser Class
	
  ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

