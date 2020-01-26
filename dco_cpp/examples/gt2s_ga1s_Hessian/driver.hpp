/*
dco/c++ v3.4.1
  - Algorithmic Differentiation by Operator Overloading in C++

COPYRIGHT
The Numerical Algorithms Group Limited and
Software and Tools for Computational Engineering @ RWTH Aachen University, 2017

This file is part of dco/c++.


*/

void t2s_a1s_Hessian(
    const vector<double>& xv,
    double& yv, 	   
    vector<vector<double> >& hessian 
);

void a2s_a1s_Hessian(
    const vector<double>& xv,
    double& yv, 	   
    vector<vector<double> >& hessian 
);

