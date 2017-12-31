#ifndef BASIS_H
#define BASIS_H

/**
 * file basis.h
 */

#include "poly.h"
#include <armadillo>
#include <cmath>
#include <time.h>
#include <iostream>
class Basis{
 public:
  int mMax;
  double br;
  double bz;
  arma::ivec nMax;

  arma::imat n_zMax;

  Basis(double,double,int,double);
  arma::vec rPart(arma::vec,int,int);
  arma::vec zPart(arma::vec,int);
  arma::mat basisFunc(int,int,int,arma::vec,arma::vec);
  arma::mat solutionref(arma::vec,arma::vec,int,int);
  arma::mat solution1(arma::vec,arma::vec,int,int);
  arma::mat solution2(arma::vec,arma::vec,int,int);

};

#endif
