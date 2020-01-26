/*
dco/c++ v3.4.1
  - Algorithmic Differentiation by Operator Overloading in C++

COPYRIGHT
The Numerical Algorithms Group Limited and
Software and Tools for Computational Engineering @ RWTH Aachen University, 2017

This file is part of dco/c++.


*/

#include <cassert>
#include <iostream>
#include <fstream>
#include <cfloat>
using namespace std;

#include "../config.h"
#include "../ad_data_view.hpp"
#include "../primal/mc.hpp"

int main(int argc, char* argv[]) {
  if(argc!=2) { cerr << "Please specify input scenario file, e.g. scenario_1.in\n"; exit(-1); }
  cout.precision(15);
  ifstream vols(argv[1]); 
  if(vols.fail()) { cerr << "Cannot open scenario file '" << argv[1] << "'\n"; exit(-1); }
  ACTIVE_INPUTS<double> X(S0,r,K,T,vols); vols.close();
  PASSIVE_INPUTS XP(N,M);
  ACTIVE_OUTPUTS<double> Y;
  ACTIVE_OUTPUTS<double> Yp;
  PASSIVE_OUTPUTS YP;

  // XM[xmsz] is an array pointing to all active inputs
  unsigned int xmsz = 4+X.sigmaSq.a.size()+X.sigmaSq.b.size();
  double **XM=new double*[xmsz];
  XM[0]=&X.S0; XM[1]=&X.r; XM[2]=&X.K; XM[3]=&X.T;
  for (unsigned int i=0;i<X.sigmaSq.a.size();i++)
    XM[i+4]=&X.sigmaSq.a[i];
  for (unsigned int i=0;i<X.sigmaSq.b.size();i++)
    XM[i+4+X.sigmaSq.a.size()]=&X.sigmaSq.b[i];

  price(X,XP,Y,YP);
  cout << "Y=" << Y.V << " in (" << Y.V-3*YP.ci << ", " << Y.V+3*YP.ci << ")" << endl;

  for (unsigned int i=0;i<xmsz;i++) {
    double p=(*XM[i]==0) ? sqrt(DBL_EPSILON) : sqrt(DBL_EPSILON)*abs(*XM[i]);
    *XM[i]+=p;
    XP.reseed(); price(X,XP,Y,YP);
    *XM[i]-=2*p;
    XP.reseed(); price(X,XP,Yp,YP);
    *XM[i]+=p;
    cout << "dY/dX[" << i << "]=" << (Y.V-Yp.V)/(2*p) << endl;
  }

  delete [] XM;
  return 0;
}

