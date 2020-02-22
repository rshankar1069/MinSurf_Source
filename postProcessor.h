/*  Simulation Sciences Laboratory
 *  WS 2019/20
 *  Chenfei Fan, Praveen Mishra, Sankrarasubramanian Ragunathan, Philipp Schleich
 *  Project 1 - "Minimal Surfaces"
 *
 *  Post-processing class -- header file
 */

#ifndef POSTPROCESS_H_INCLUDED
#define POSTPROCESS_H_INCLUDED

#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkStructuredGrid.h>
#include <vtkSmartPointer.h>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include "inputParser.h"
#include "cartesianGrid.h"
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <cstdlib>
#include <sys/types.h>
#include <sys/stat.h>

#define _USE_MATH_DEFINES

// Function to write the structured VTK file for visualization
template<class mType,class dType>
void structuredGridWriter(int iterationIndex,mType z);

// Function to write the residual
template <class dType>
void residualWriter(int iterationIndex, dType res);

// Function to compare Analytical solution with Numerical solution
template <class mType, class dType>
void compareAnalyticSoln(int opt, mType z);

// Function to calculate the l2-Euclidean Norm
template <class mType, class dType>
dType l2Euclidean(int N, mType z, std::vector<dType> zAnalytical);

//Function to calculate the Max Norm
template <class mType, class dType>
dType maxNorm(int N, mType z, std::vector<dType> zAnalytical);

#endif //POSTPROCESS_H_INCLUDED
