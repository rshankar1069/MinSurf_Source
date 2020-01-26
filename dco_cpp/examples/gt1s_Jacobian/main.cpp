/*
dco/c++ v3.4.1
  - Algorithmic Differentiation by Operator Overloading in C++

COPYRIGHT
The Numerical Algorithms Group Limited and
Software and Tools for Computational Engineering @ RWTH Aachen University, 2017

This file is part of dco/c++.


*/

#include<iostream>
#include<vector>
using namespace std;

#include "driver.hpp"

int main() {
  cout.precision(15);
  const int m=2, n=4; 
  vector<double> x(n), y(m);
  vector<vector<double> > J(m,x); 
  for (int i=0;i<n;i++) x[i]=1; 
  driver(x,y,J);
  for (int i=0;i<m;i++) 
    cout << "y[" << i << "]=" << y[i] << endl;
  for (int j=0;j<m;j++) 
    for (int i=0;i<n;i++) 
      cout << "F'[" << j << "][" << i << "]=" << J[j][i] << endl;
  return 0;
}

