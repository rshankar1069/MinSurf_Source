/*
dco/c++ v3.4.1
  - Algorithmic Differentiation by Operator Overloading in C++

COPYRIGHT
The Numerical Algorithms Group Limited and
Software and Tools for Computational Engineering @ RWTH Aachen University, 2017

This file is part of dco/c++.


*/

#include <vector>
using namespace std;

#include "dco.hpp" 
using namespace dco;

typedef long double DCO_VALUE_TYPE;
typedef double DCO_PARTIAL_TYPE;
typedef float DCO_ADJOINT_TYPE;
typedef ga1s<DCO_VALUE_TYPE,DCO_PARTIAL_TYPE,DCO_ADJOINT_TYPE> DCO_MODE;
typedef DCO_MODE::type DCO_TYPE;
typedef DCO_MODE::tape_t DCO_TAPE_TYPE;

#include "f.hpp"
		   
void driver(
    const vector<double>& xv, vector<double>& xa,     
    vector<double>& yv, vector<double>& ya      
) {
  DCO_MODE::global_tape=DCO_TAPE_TYPE::create();
  size_t n=xv.size(), m=yv.size();
  vector<DCO_TYPE> x(n), y(m);
  for (size_t i=0;i<n;i++) {
    x[i]=xv[i];
    DCO_MODE::global_tape->register_variable(x[i]);
  }
  f(x,y);
  for (size_t i=0;i<m;i++) {
    DCO_MODE::global_tape->register_output_variable(y[i]);
    yv[i]=static_cast<float>(value(y[i])); derivative(y[i])=static_cast<float>(ya[i]);
  }
  for (size_t i=0;i<n;i++) derivative(x[i])=static_cast<float>(xa[i]);

  DCO_MODE::global_tape->write_to_dot();
  DCO_MODE::global_tape->interpret_adjoint();
  for (size_t i=0;i<n;i++) xa[i]=derivative(x[i]);
  for (size_t i=0;i<m;i++) ya[i]=derivative(y[i]);
  DCO_TAPE_TYPE::remove(DCO_MODE::global_tape);
}

