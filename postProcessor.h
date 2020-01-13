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

#define _USE_MATH_DEFINES

// Function to write the structured VTK file of the final solution for visualization
template<class dType>
void structuredGridWriter(std::string filename,int iterationIndex);

#endif //POSTPROCESS_H_INCLUDED
