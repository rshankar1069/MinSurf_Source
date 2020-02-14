#include"cartesianGrid.h"
#include"pch.h"
#include "iostream.h"
#include "cmath.h"

// -------------------------------------------------------------------------------------------------
template <class dType, class listType>
void cartesianGrid<dType, listType>::setupGrid() {

	setNoGridPoints(101);
	setGridSpacing();

	setBdryNodes();
	setInnerNodes();

}



//Lexicographical Ordering -------------------------------------------------------------------------------------------------

template <class dType, class listType>
void cartesianGrid<dType, listType>::setBdryNodes() {

	// Set first and last row
	for (int i = 0; i < noGridPoints; i++) {
		// First row
		bdryNodeList.bottom.push_back(i);
		// Last row
		bdryNodeList.top.push_back(noGridPoints * (noGridPoints - 1) + i);
		//         std::cout << "me: " << i << std::endl;
	}

	// outer edges
	for (int i = noGridPoints; i < noGridPoints * (noGridPoints - 1); i += noGridPoints) {
		bdryNodeList.left.push_back(i);
		bdryNodeList.right.push_back(i + noGridPoints - 1);
	}

	std::sort(bdryNodeList.bottom.begin(), bdryNodeList.bottom.end());
	std::sort(bdryNodeList.right.begin(), bdryNodeList.right.end());
	std::sort(bdryNodeList.top.begin(), bdryNodeList.top.end());
	std::sort(bdryNodeList.left.begin(), bdryNodeList.left.end());


}


// -------------------------------------------------------------------------------------------------
// Inner nodes by lexicographical ordering
template <class dType, class listType>
void cartesianGrid<dType, listType>::setInnerNodes() {

	// Set inner nodes (assuming structured cartesian grid)
	for (int i = 1; i < noGridPoints - 1; i++) {
		for (int j = 1; j < noGridPoints - 1; j++) {
			innerNodeList.push_back(i + j * noGridPoints);
		}
	}
	// Sort list (works for vector)
	std::sort(innerNodeList.begin(), innerNodeList.end());


}

// Get x coordinate of a certain node (index)
template <class dType, class listType>
dType cartesianGrid<dType, listType>::getCooX(int index) {
	int column = index % noGridPoints;
	dType x = ((dType)column) / (noGridPoints - 1);
	return x;
}


// Get y coordinate of a certain node (index)
template <class dType, class listType>
dType cartesianGrid<dType, listType>::getCooY(int index) {
	int row = index / noGridPoints;
	dType y = ((dType)row) / (noGridPoints - 1);
	return y;
}
}