/*
dco/c++ v3.4.1
  - Algorithmic Differentiation by Operator Overloading in C++

COPYRIGHT
The Numerical Algorithms Group Limited and
Software and Tools for Computational Engineering @ RWTH Aachen University, 2017

This file is part of dco/c++.


*/

#include <cassert>
#include <cstdlib>
#include <iostream>
using namespace std;

#include "driver.hpp"

int main() {
    int n=10; cout.precision(15);
    double xv=2.1,xa=1;
    double pv=1,pa=0;
    driver(n,xv,xa,pv,pa);
    cout << "x=" << xv << endl;
    cout << "x_{(1)}=" << xa << endl;
    cout << "p_{(1)}=" << pa << endl;
    return 0;
}


