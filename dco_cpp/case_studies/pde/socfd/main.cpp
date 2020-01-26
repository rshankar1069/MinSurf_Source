/*
dco/c++ v3.4.1
  - Algorithmic Differentiation by Operator Overloading in C++

COPYRIGHT
The Numerical Algorithms Group Limited and
Software and Tools for Computational Engineering @ RWTH Aachen University, 2017

This file is part of dco/c++.


*/

#include <assert.h>
#include <iostream>
#include <fstream>
#include <cfloat>
using namespace std;

#include "../config.h"
#include "../ad_data_view.hpp"
#include "../primal/pdesolver.hpp"

int main(int argc, char* argv[]) {
  if(argc!=2) { cerr << "Please specify input scenario file, e.g. scenario_1.in\n"; exit(-1); }
  cout.precision(15);
  ifstream vols(argv[1]);
  if(vols.fail()) { cerr << "Cannot open scenario file '" << argv[1] << "'\n"; exit(-1); }
  ACTIVE_INPUTS<double> X(S0,r,K,T,vols);
  vols.close();
  PASSIVE_INPUTS XP(N,M);

  ACTIVE_OUTPUTS<double> Y;
  ACTIVE_OUTPUTS<double> Yp;

  // XM[xmsz] is an array pointing to all active inputs
  unsigned int xmsz = 4+X.sigmaSq.a.size()+X.sigmaSq.b.size();
  double **XM=new double*[xmsz];
  XM[0]=&X.S0; XM[1]=&X.r; XM[2]=&X.K; XM[3]=&X.T;
  for (unsigned int i=0;i<X.sigmaSq.a.size();i++)
    XM[i+4]=&X.sigmaSq.a[i];
  for (unsigned int i=0;i<X.sigmaSq.b.size();i++)
    XM[i+4+X.sigmaSq.a.size()]=&X.sigmaSq.b[i];

  price(X,XP,Y);
  cout << "Y=" << Y.V << endl;

  for (unsigned int i=0;i<xmsz;i++) {
    double pi=(*XM[i]==0) ? sqrt(sqrt(DBL_EPSILON)) : sqrt(sqrt(DBL_EPSILON))*abs(*XM[i]);
    for (unsigned int j=0;j<=i;j++) {
      double pj=(*XM[j]==0) ? sqrt(sqrt(DBL_EPSILON)) : sqrt(sqrt(DBL_EPSILON))*abs(*XM[j]);
      *XM[i]+=pi; *XM[j]+=pj;
      price(X,XP,Y);
      *XM[i]-=2*pi; 
      price(X,XP,Yp); 
      Y.V-=Yp.V;
      *XM[i]+=2*pi; *XM[j]-=2*pj;
      price(X,XP,Yp); 
      Y.V-=Yp.V;
      *XM[i]-=2*pi; 
      price(X,XP,Yp); 
      Y.V+=Yp.V;
      cout << "d^2Y/dX[" << i << "]dX[" << j << "]=" 
	   << Y.V/(4*pi*pj) << endl;
      *XM[i]+=pi; *XM[j]+=pj;
    }
  }

  delete [] XM;
  return 0;
}

