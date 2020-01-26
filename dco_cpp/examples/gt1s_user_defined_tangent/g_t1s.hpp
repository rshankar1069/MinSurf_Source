/*
dco/c++ v3.4.1
  - Algorithmic Differentiation by Operator Overloading in C++

COPYRIGHT
The Numerical Algorithms Group Limited and
Software and Tools for Computational Engineering @ RWTH Aachen University, 2017

This file is part of dco/c++.


*/

template<>
void g<gt1s<double>::type>(
    const gt1s<double>::type& p, 
    gt1s<double>::type& x
) {
  double pv,xv,pt1,xt1;
  xv=value(x);
  pv=value(p);
  g(pv,xv);
  pt1=derivative(p);
  xt1=pt1/(2*xv);
  value(x)=xv;
  derivative(x)=xt1;
}

