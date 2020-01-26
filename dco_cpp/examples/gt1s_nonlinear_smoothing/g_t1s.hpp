/*
dco/c++ v3.4.1
  - Algorithmic Differentiation by Operator Overloading in C++

COPYRIGHT
The Numerical Algorithms Group Limited and
Software and Tools for Computational Engineering @ RWTH Aachen University, 2017

This file is part of dco/c++.


*/

#include<iostream>

template<>
void g<dco::gt1s<double>::type>(dco::gt1s<double>::type& x) {
	double xv = dco::value(x);
  x=-(.5+atan(x)/atan(1.));
  std::cout << x << std::endl;
  if (xv>0)
    xv=0;
  else 
    xv=1;
  dco::value(x) = xv;
}

